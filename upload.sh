#!/bin/bash
git add .
git commit -m "${1:-Update blog}"
git push origin main