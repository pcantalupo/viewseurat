# Ship Command
Merge current branch to main and clean up.

Confirm the branch name before starting, then launch a Bash subagent Task to do all of the following steps (this allows the user to Ctrl+B to background it):

* Run `Rscript -e 'devtools::document()'` to update NAMESPACE and .Rd files
* Check for uncommitted changes - if any, commit them with a concise message
* Switch to main and pull latest
* Merge the feature branch into main (fast-forward if possible)
* Push main to origin
* Delete the local feature branch
* Delete the remote feature branch (origin)
* Prune any stale remote-tracking branches

If there are merge conflicts, stop and report.

