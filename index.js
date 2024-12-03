function showModal() {
    document.getElementById('cvModal').style.display = 'block';
    document.getElementById('cvModalOverlay').style.display = 'block';
}

function hideModal() {
    document.getElementById('cvModal').style.display = 'none';
    document.getElementById('cvModalOverlay').style.display = 'none';
}

function playSound() {
    var audio = document.getElementById('sound');
    audio.play();
}
