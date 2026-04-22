#!/bin/bash
set -euo pipefail

# fail if repo is dirty:
git diff --quiet && git diff --cached --quiet || {
  echo "Repository has uncommitted changes"
  exit 1
}

# fail on untracked files, also fail on ignored to reduce risk of missing edits for a release
[ -z "$(git status --porcelain --ignored)" ] || {
  echo "Repository has changes (including ignored files)"
  exit 1
}

echo merging to releases

# integration
git fetch origin releases:releases

TMP="${TMP:-/tmp/}"

# check version follows the last release
version="$(cat version.txt)"
last_release_tmpfile=$(mktemp $TMP/oncoscanner_last_release_tmpfile.XXXXXX)
git show releases:version.txt > $last_release_tmpfile #../last_release_version.txt
python scripts/version_tool.py follows $last_release_tmpfile version.txt #../last_release_version.txt version.txt
rm $last_release_tmpfile

# instead full clone and then
# ensure tip of main branch
git fetch

if [ "$(git rev-parse HEAD)" != "$(git rev-parse main)" ]; then
	echo "Not at tip of branch (ahead or behind upstream)"
	exit 1
fi

echo checking out releases branch
git checkout releases
mkdir -p metadata

GIT_EDITOR=true git merge --no-ff main --into releases --no-commit || true

# save last commit to file (as version metadata)
git log -1 --pretty=format:"%H" main > metadata/prev_commit_hash.txt
git describe --always --abbrev=999 > metadata/git_describe.txt
git add metadata/prev_commit_hash.txt metadata/git_describe.txt

version_line="manifest.version=$version"
sed -i '' "s/^manifest\.version.*/$version_line/" modules/oncoscanner/user.config

git commit -m 'create release'

# will conflict on version.txt
# create annotated tag for the release
git tag -a "$version" -m "$version"

#git checkout main version.txt
#GIT_EDITOR=true git merge --continue
#git push

echo checking out back:
git checkout main

echo release has been "made" successfully,
echo now you can git push and git push --tags
