---
title: "SHAN: Object-Level Privacy Detection via Inference on Scene Heterogeneous Graph"
authors:
- admin
- Robert Ford
author_notes:
- "Equal contribution"
- "Equal contribution"
date: "2024-03-01T00:00:00Z"
doi: ""

# Schedule page publish date (NOT publication's date).
publishDate: "2024-03-01T00:00:00Z"

# Publication type.
# Accepts a single type but formatted as a YAML list (for Hugo requirements).
# Enter a publication type from the CSL standard.
publication_types: ["article-journal"]

# Publication name and optional abbreviated publication name.
publication: "*Journal of Source Themes, 1*(1)"
publication_short: ""

abstract: With the rise of social platforms, protecting privacy has become an important issue. Privacy object detection aims to accurately locate private objects in images. It is the foundation of safeguarding individuals' privacy rights and ensuring responsible data handling practices in the digital age. Since privacy of object is not shift-invariant, the essence of the privacy object detection task is inferring object privacy based on scene information. However, privacy object detection has long been studied as a subproblem of common object detection tasks. Therefore, existing methods suffer from serious deficiencies in accuracy, generalization, and interpretability. Moreover, creating large-scale privacy datasets is difficult due to legal constraints and existing privacy datasets lack label granularity. The granularity of existing privacy detection methods remains limited to the image level. To address the above two issues, we introduce two benchmark datasets for object-level privacy detection and propose SHAN, Scene Heterogeneous graph Attention Network, a model constructs a scene heterogeneous graph from an image and utilizes self-attention mechanisms for scene inference to obtain object privacy. Through experiments, we demonstrated that SHAN performs excellently in privacy object detection tasks, with all metrics surpassing those of the baseline model.
# Summary. An optional shortened abstract.
summary: With the rise of social platforms, protecting privacy has become an important issue. Privacy object detection aims to accurately locate private objects in images. It is the foundation of safeguarding individuals' privacy rights and ensuring responsible data handling practices in the digital age. Since privacy of object is not shift-invariant, the essence of the privacy object detection task is inferring object privacy based on scene information. However, privacy object detection has long been studied as a subproblem of common object detection tasks. Therefore, existing methods suffer from serious deficiencies in accuracy, generalization, and interpretability. Moreover, creating large-scale privacy datasets is difficult due to legal constraints and existing privacy datasets lack label granularity. The granularity of existing privacy detection methods remains limited to the image level. To address the above two issues, we introduce two benchmark datasets for object-level privacy detection and propose SHAN, Scene Heterogeneous graph Attention Network, a model constructs a scene heterogeneous graph from an image and utilizes self-attention mechanisms for scene inference to obtain object privacy. Through experiments, we demonstrated that SHAN performs excellently in privacy object detection tasks, with all metrics surpassing those of the baseline model.
tags:
- Source Themes
featured: false

# links:
# - name: ""
#   url: ""
url_pdf: http://arxiv.org/pdf/1512.04133v1
url_code: 'https://github.com/HugoBlox/hugo-blox-builder'
url_dataset: ''
url_poster: ''
url_project: ''
url_slides: ''
url_source: ''
url_video: ''

# Featured image
# To use, add an image named `featured.jpg/png` to your page's folder. 
image:
  caption: 'Image credit: [**Unsplash**](https://unsplash.com/photos/jdD8gXaTZsc)'
  focal_point: ""
  preview_only: false

# Associated Projects (optional).
#   Associate this publication with one or more of your projects.
#   Simply enter your project's folder or file name without extension.
#   E.g. `internal-project` references `content/project/internal-project/index.md`.
#   Otherwise, set `projects: []`.
projects: []

# Slides (optional).
#   Associate this publication with Markdown slides.
#   Simply enter your slide deck's filename without extension.
#   E.g. `slides: "example"` references `content/slides/example/index.md`.
#   Otherwise, set `slides: ""`.
slides: example
---
