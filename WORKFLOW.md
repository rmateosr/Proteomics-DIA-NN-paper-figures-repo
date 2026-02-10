# Git Workflow Guide

**Repository:** https://github.com/rmateosr/Proteomics-DIA-NN-paper-figures-repo
**Local directory:** `C:\Users\Raul\Dropbox\Papers\DIA-NN\Figures\scripts`

---

## How This Repository Was Created

The following steps were performed on 2026-02-10 using the GitHub CLI (`gh`) from WSL2:

### 1. Authenticate with GitHub
```bash
gh auth login --web -p https
```
This opened a browser for OAuth authentication and logged in as `rmateosr`.

### 2. Initialize a git repository in the scripts folder
```bash
cd "C:\Users\Raul\Dropbox\Papers\DIA-NN\Figures\scripts"
git init
git branch -m main
```

### 3. Configure git identity
```bash
git config user.name "rmateosr"
git config user.email "rmateosr@users.noreply.github.com"
```

### 4. Stage and commit all R scripts
```bash
git add *.R
git commit -m "Add R scripts for DIA-NN proteomics paper figures"
```

### 5. Create the remote GitHub repository and push
```bash
gh repo create rmateosr/Proteomics-DIA-NN-paper-figures-repo \
  --public \
  --description "R scripts for generating figures in the DIA-NN proteomics paper" \
  --source . \
  --push
```

### 6. Set up credential helper for future pushes
```bash
gh auth setup-git
```

This created the public repository at https://github.com/rmateosr/Proteomics-DIA-NN-paper-figures-repo with the working scripts directory as the local repo — no separate folder needed.

---

## Setup

Your working scripts folder is already a git repository linked to GitHub. You edit scripts in place and use git to track and push changes.

---

## Daily Workflow

### 1. After editing a script

Open a terminal (Git Bash, WSL, or Command Prompt) in `C:\Users\Raul\Dropbox\Papers\DIA-NN\Figures\scripts` and run:

```bash
# See what changed
git status

# Stage the modified file(s)
git add <filename>.R

# Commit with a short description of what you changed
git commit -m "Update KRAS figure colors and labels"

# Push to GitHub
git push
```

### 2. Adding a new script

```bash
# Stage the new file
git add new_script.R

# Commit
git commit -m "Add new analysis script for XYZ"

# Push
git push
```

### 3. Checking what has changed before committing

```bash
# See which files were modified
git status

# See the actual line-by-line changes
git diff

# See changes for a specific file
git diff <filename>.R
```

### 4. Viewing commit history

```bash
# See recent commits
git log --oneline -10
```

---

## Tips

- **Commit often:** Small, frequent commits are better than one large commit. Each commit should represent one logical change (e.g., "Fix axis labels in TP53 figure").
- **Write clear messages:** Future you will thank present you. Use messages like "Add CTNNB1 panel to PDX hotspot figure" rather than "update".
- **Push regularly:** `git push` uploads your commits to GitHub. Until you push, changes only exist locally.
- **Don't worry about breaking things:** Git keeps full history. You can always go back to any previous version with `git log` and `git checkout`.

---

## Common Scenarios

### "I want to undo my last uncommitted changes to a file"
```bash
git checkout -- <filename>.R
```
This restores the file to the last committed version. Uncommitted changes will be lost.

### "I want to see an older version of a script"
```bash
# Find the commit hash from the log
git log --oneline <filename>.R

# View the file at that commit
git show <commit-hash>:<filename>.R
```

### "I committed something by mistake"
```bash
# Undo the last commit but keep the changes in your files
git reset --soft HEAD~1
```

---

## Using from RStudio

If you use RStudio, you can also manage git from its built-in Git pane:

1. Open the project in RStudio (File > Open Project, navigate to the scripts folder)
2. The **Git** tab (top-right pane) shows changed files
3. Check the boxes to stage files
4. Click **Commit**, write a message, and click **Commit**
5. Click **Push** (green up arrow) to upload to GitHub
