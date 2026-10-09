(() => {
  'use strict';

  const CONFIG = window.HUMBLE_EVAL_CONFIG || {};
  const ANNOTATOR = Number(document.body.dataset.annotator);
  const STORE_KEY = `humble-eval-a${ANNOTATOR}`;
  const LETTERS = ['A', 'B', 'C', 'D', 'E'];
  const app = document.getElementById('app');

  let items = [];
  let state = null;
  let view = 'intro';   // intro | item | done
  let flushing = false;

  // Active (visible) time on the current item.
  let activeSince = null;
  let activeAccum = 0;
  let shownWall = Date.now();   // wall-clock fallback if visibility reporting is off

  // ---------- persistence ----------

  function newState() {
    return {
      annotator: ANNOTATOR,
      client_id: Math.random().toString(36).slice(2, 10),
      started_at: null,
      finished_at: null,
      current: 0,
      responses: {},   // uid -> response record
      queue: [],       // events waiting to be sent
      final_queued: false,
    };
  }

  function load() {
    try {
      const raw = localStorage.getItem(STORE_KEY);
      return raw ? JSON.parse(raw) : null;
    } catch (e) {
      return null;
    }
  }

  function save() {
    try {
      localStorage.setItem(STORE_KEY, JSON.stringify(state));
    } catch (e) { /* storage unavailable: server sync still works */ }
  }

  // ---------- server sync ----------

  function enqueue(event) {
    if (!CONFIG.endpoint) return;
    state.queue.push({
      ...event,
      event_id: `${state.client_id}-${Date.now()}-${Math.random().toString(36).slice(2, 6)}`,
      annotator: ANNOTATOR,
      client_id: state.client_id,
    });
    save();
    flush();
  }

  async function flush() {
    if (!CONFIG.endpoint || flushing) return;
    flushing = true;
    try {
      while (state.queue.length) {
        const ev = state.queue[0];
        const res = await fetch(CONFIG.endpoint, {
          method: 'POST',
          headers: { 'Content-Type': 'text/plain;charset=utf-8' },
          body: JSON.stringify({ ...ev, token: CONFIG.token }),
        });
        const text = await res.text();
        if (!text.startsWith('ok')) throw new Error(text);
        state.queue.shift();
        save();
        renderStatus();
      }
    } catch (e) {
      // Network or server problem: keep the queue and retry later.
    } finally {
      flushing = false;
      renderStatus();
    }
  }

  setInterval(flush, 15000);
  window.addEventListener('online', flush);

  // ---------- timing ----------

  function startTimer() {
    activeAccum = 0;
    activeSince = document.hidden ? null : performance.now();
    shownWall = Date.now();
  }

  function elapsed() {
    return activeAccum + (activeSince === null ? 0 : performance.now() - activeSince);
  }

  document.addEventListener('visibilitychange', () => {
    if (view !== 'item') return;
    if (document.hidden && activeSince !== null) {
      activeAccum += performance.now() - activeSince;
      activeSince = null;
    } else if (!document.hidden && activeSince === null) {
      activeSince = performance.now();
    }
  });

  // ---------- actions ----------

  function answeredCount() {
    return Object.keys(state.responses).length;
  }

  function firstUnanswered() {
    return items.findIndex(it => !state.responses[it.uid]);
  }

  function answer(letter) {
    const it = items[state.current];
    const now = new Date().toISOString();
    const spent = Math.round(elapsed());
    const wall = Date.now() - shownWall;
    const prev = state.responses[it.uid];

    if (prev && prev.answer === letter) {
      goNext();
      return;
    }

    const rec = prev ? {
      ...prev,
      answer: letter,
      n_changes: prev.n_changes + 1,
      total_time_ms: prev.total_time_ms + spent,
      total_wall_ms: (prev.total_wall_ms || 0) + wall,
      answered_at: now,
    } : {
      index: it.index,
      uid: it.uid,
      answer: letter,
      first_time_ms: spent,
      total_time_ms: spent,
      first_wall_ms: wall,
      total_wall_ms: wall,
      n_changes: 0,
      first_answered_at: now,
      answered_at: now,
    };
    state.responses[it.uid] = rec;
    enqueue({ type: 'answer', ...rec });
    render();
    setTimeout(goNext, 180);
  }

  function goTo(i) {
    state.current = Math.max(0, Math.min(items.length - 1, i));
    save();
    view = 'item';
    render();
  }

  function goNext() {
    if (state.current < items.length - 1) {
      goTo(state.current + 1);
    } else if (firstUnanswered() >= 0) {
      goTo(firstUnanswered());
    } else {
      finish();
    }
  }

  function finish() {
    if (firstUnanswered() >= 0) {
      goTo(firstUnanswered());
      return;
    }
    if (!state.final_queued) {
      state.finished_at = new Date().toISOString();
      state.final_queued = true;
      enqueue({
        type: 'final',
        n_answered: answeredCount(),
        started_at: state.started_at,
        finished_at: state.finished_at,
        user_agent: navigator.userAgent,
        responses: items.map(it => state.responses[it.uid]),
      });
    }
    view = 'done';
    render();
  }

  function downloadBackup() {
    const blob = new Blob([JSON.stringify({
      annotator: ANNOTATOR,
      client_id: state.client_id,
      started_at: state.started_at,
      finished_at: state.finished_at,
      downloaded_at: new Date().toISOString(),
      n_answered: answeredCount(),
      user_agent: navigator.userAgent,
      responses: items.map(it => state.responses[it.uid] || null),
    }, null, 2)], { type: 'application/json' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `annotator_${ANNOTATOR}_results.json`;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  // ---------- rendering ----------

  function esc(s) {
    return String(s).replace(/[&<>"']/g, c => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
    }[c]));
  }

  function renderStatus() {
    const n = answeredCount();
    const prog = document.getElementById('progress');
    const bar = document.getElementById('bar');
    const sync = document.getElementById('sync');
    if (prog) prog.textContent = `${n} / ${items.length} answered`;
    if (bar) bar.style.width = `${items.length ? (100 * n / items.length) : 0}%`;
    updateDoneSync();
    if (!sync) return;
    if (!CONFIG.endpoint) {
      sync.textContent = n ? 'Progress saved in this browser' : '';
      sync.className = 'sync';
    } else if (state.queue.length === 0) {
      sync.textContent = n ? 'All answers saved ✓' : '';
      sync.className = 'sync ok';
    } else {
      sync.textContent = `Saving… (${state.queue.length} pending)`;
      sync.className = 'sync warn';
    }
  }

  function renderIntro() {
    const resuming = answeredCount() > 0;
    app.innerHTML = `
      <div class="card intro">
        <h2>Image question study</h2>
        <p>You will see <b>${items.length}</b> multiple-choice questions about images.</p>
        <ul>
          <li>Look at the image carefully and choose the option that is correct for the image.</li>
          <li>If none of the options A–D is correct, choose <b>E (None of the above)</b>.</li>
          <li>Please answer based on what you can see in the image. There is no time limit.</li>
          <li>Your progress is saved automatically. You can close the page and continue later in the same browser.</li>
          <li>When you finish, click <b>Download results</b> and send the file to the organizer.</li>
          <li>Keyboard shortcuts: <b>A–E</b> or <b>1–5</b> to answer, <b>←</b>/<b>→</b> to move.</li>
        </ul>
        <button class="btn primary" id="start">${resuming ? `Continue (${answeredCount()} / ${items.length} done)` : 'Start'}</button>
      </div>`;
    document.getElementById('start').onclick = () => {
      if (!state.started_at) {
        state.started_at = new Date().toISOString();
        save();
      }
      if (answeredCount() === items.length) finish();
      else goTo(resuming && firstUnanswered() >= 0 ? firstUnanswered() : state.current);
    };
  }

  function renderItem() {
    const it = items[state.current];
    const chosen = state.responses[it.uid]?.answer;
    const options = LETTERS.filter(l => l in it.options).map(l => `
      <button class="option${chosen === l ? ' selected' : ''}" data-letter="${l}">
        <span class="letter">${l}</span><span>${esc(it.options[l])}</span>
      </button>`).join('');
    const allDone = firstUnanswered() < 0;

    app.innerHTML = `
      <div class="card item">
        <div class="image-wrap"><img src="${esc(it.image)}" alt="Question image"></div>
        <div>
          <div class="qnum">Question ${state.current + 1} of ${items.length}</div>
          <p class="question">${esc(it.question)}</p>
          <div class="options">${options}</div>
          <div class="nav">
            <button class="btn" id="prev" ${state.current === 0 ? 'disabled' : ''}>← Previous</button>
            <button class="btn" id="next" ${state.current === items.length - 1 ? 'disabled' : ''}>Next →</button>
            ${allDone ? '<button class="btn primary" id="finish">Finish</button>' : ''}
          </div>
          <div class="hint">Keys: A–E or 1–5 to answer · ← / → to move</div>
        </div>
      </div>`;

    app.querySelectorAll('.option').forEach(b => {
      b.onclick = () => answer(b.dataset.letter);
    });
    document.getElementById('prev').onclick = () => goTo(state.current - 1);
    document.getElementById('next').onclick = () => goTo(state.current + 1);
    const fin = document.getElementById('finish');
    if (fin) fin.onclick = finish;

    // Preload the next image.
    const next = items[state.current + 1];
    if (next) new Image().src = next.image;
    startTimer();
  }

  function renderDone() {
    app.innerHTML = `
      <div class="card done">
        <p class="big">All ${items.length} questions answered. Thank you!</p>
        <p id="done-sync"></p>
        <p>${CONFIG.endpoint
          ? '<button class="small-link" id="backup">Download a backup copy</button>'
          : '<button class="btn primary" id="backup">Download results</button>'}</p>
        <p><button class="btn" id="review">Review my answers</button></p>
      </div>`;
    document.getElementById('review').onclick = () => goTo(0);
    document.getElementById('backup').onclick = downloadBackup;
    updateDoneSync();
  }

  function updateDoneSync() {
    const el = document.getElementById('done-sync');
    if (!el) return;
    if (!CONFIG.endpoint) {
      el.textContent = `Please download your results and send the file (annotator_${ANNOTATOR}_results.json) to the organizer.`;
    } else if (state.queue.length) {
      el.textContent = 'Saving your answers… please keep this page open for a moment.';
    } else {
      el.textContent = 'Your answers have been saved. You can close this page.';
    }
  }

  function render() {
    if (view === 'intro') renderIntro();
    else if (view === 'item') renderItem();
    else renderDone();
    renderStatus();
  }

  document.addEventListener('keydown', e => {
    if (view !== 'item' || e.metaKey || e.ctrlKey || e.altKey) return;
    const k = e.key.toUpperCase();
    const idx = LETTERS.indexOf(k) >= 0 ? LETTERS.indexOf(k) : '12345'.indexOf(k);
    if (idx >= 0 && LETTERS[idx] in items[state.current].options) {
      answer(LETTERS[idx]);
    } else if (e.key === 'ArrowLeft') {
      goTo(state.current - 1);
    } else if (e.key === 'ArrowRight') {
      goTo(state.current + 1);
    }
  });

  // ---------- boot ----------

  fetch(`data/annotator_${ANNOTATOR}.json`)
    .then(r => r.json())
    .then(data => {
      items = data.items;
      state = load() || newState();
      save();
      flush();
      render();
    })
    .catch(() => {
      app.innerHTML = '<div class="card intro"><p>Could not load the questions. Please refresh the page.</p></div>';
    });
})();
