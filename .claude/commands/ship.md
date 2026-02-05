Merge current branch to main and clean up:

1. Check for uncommitted changes - if any, commit them with a concise message
2. Switch to main and pull latest
3. Merge the feature branch into main (fast-forward if possible)
4. Push main to origin
5. Delete the local feature branch
6. Delete the remote feature branch (origin)
7. Prune any stale remote-tracking branches

Confirm the branch name before starting. If there are merge conflicts, stop and report.

