#!/usr/bin/env bash
# Publish the tracked wiki/ files (pages and wiki/assets/ images) to the repository's GitHub Wiki after
# docs_check.sh accepts them. The wiki must already exist: create its first page in GitHub once. Files the wiki has
# and wiki/ no longer tracks are removed one named file at a time (git rm), never by a recursive delete; the clone
# is left in RUNNER_TEMP.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

repository=${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is required}
server_url=${GITHUB_SERVER_URL:-https://github.com}
: "${WIKI_TOKEN:?WIKI_TOKEN is required}"
source_sha=$(git rev-parse HEAD)
bash .github/scripts/docs_check.sh

work=$(mktemp -d "${RUNNER_TEMP:-/var/tmp}/z-fastq-wiki.XXXXXX")
git ls-files -z -- wiki | sed -z 's,^wiki/,,' > "$work/files"

export GIT_CONFIG_COUNT=1
export GIT_CONFIG_KEY_0="http.${server_url%/}/.extraheader"
GIT_CONFIG_VALUE_0="AUTHORIZATION: basic $(printf 'x-access-token:%s' "$WIKI_TOKEN" | base64 | tr -d '\n')"
export GIT_CONFIG_VALUE_0
export GIT_TERMINAL_PROMPT=0

wiki_repo="$work/wiki"
if ! git clone --quiet --depth 1 "${server_url%/}/${repository}.wiki.git" "$wiki_repo"; then
    echo "cannot clone the GitHub Wiki; enable it and create its first page in GitHub" >&2
    exit 1
fi
wiki_branch=$(git -C "$wiki_repo" symbolic-ref --quiet --short HEAD)

while IFS= read -r -d '' name; do
    grep -qzFx -- "$name" "$work/files" || git -C "$wiki_repo" rm --quiet -- "$name"
done < <(git -C "$wiki_repo" ls-files -z)
while IFS= read -r -d '' name; do
    mkdir -p "$wiki_repo/$(dirname -- "$name")"
    cp -- "wiki/$name" "$wiki_repo/$name"
done < "$work/files"

git -C "$wiki_repo" add --all
if git -C "$wiki_repo" diff --cached --quiet --exit-code; then
    echo "wiki is already synchronized"
    exit 0
fi
git -C "$wiki_repo" diff --cached --name-status
git -C "$wiki_repo" config user.name 'github-actions[bot]'
git -C "$wiki_repo" config user.email '41898282+github-actions[bot]@users.noreply.github.com'
git -C "$wiki_repo" commit --quiet -m "docs: sync wiki from ${source_sha:0:12}"
git -C "$wiki_repo" push --quiet origin "HEAD:$wiki_branch"

local_sha=$(git -C "$wiki_repo" rev-parse HEAD)
remote_sha=$(git -C "$wiki_repo" ls-remote origin "refs/heads/$wiki_branch" | awk 'NR == 1 { print $1 }')
if [[ "$remote_sha" != "$local_sha" ]]; then
    echo "wiki push verification failed: local=$local_sha remote=${remote_sha:-missing}" >&2
    exit 1
fi
echo "published wiki commit $local_sha from source $source_sha"
