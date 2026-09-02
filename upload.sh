#!/bin/bash
git add .
git commit -m "Auto upload: $(date +'%Y-%m-%d %H:%M:%S')"
git push origin main