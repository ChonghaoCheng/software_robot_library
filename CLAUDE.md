# software_robot_library

## ARA: agent-native research artifacts

This project records its research in an ARA artifact
(https://github.com/ARA-Labs/Agent-Native-Research-Artifact).

**This project has two artifact locations. They are not equivalent:**

- `/home/eric/Workspace/MPCC_ARA/` — **canonical**. Its own git repo, carrying the
  full `evidence/` + `logic/` + `trace/` history. All new records go here.
- `ara/` (inside this repo) — a partial mirror, currently **stale**: it stops at
  claim C07 and is missing the F2-01/F2-02 evidence, and 7 shared files have
  diverged. Treat it as read-only. Do not write research records here, and do not
  reconcile the two without asking.

Wherever a skill below takes `<ara-dir>`, that means `/home/eric/Workspace/MPCC_ARA/`.

Route work to the matching ARA skill — invoke these yourself, without being asked:

- `/research-manager` — trigger whenever a research milestone lands: an
  experiment finishes, a decision is made, a hypothesis is confirmed or killed,
  a dead end is hit, a direction pivots, user's input. This holds equally in autonomous runs
  (loops, heartbeats, long experiments) where the user gives no input at all —
  crystallize the insight at the milestone. It
  records what just happened (decisions, experiments, dead ends, claims) into
  the artifact. Skip when nothing research-significant happened (greetings, pure formatting).
- `/compiler <path>` — when turning an existing paper, repo, logs, or notes into
  a structured artifact.
- `/rigor-reviewer <dir>` — before trusting, publishing, or submitting an artifact.
- `/research-visualizer <ara-dir>` — to inspect the research trajectory as an
  interactive process map (add `--serve` for a live local viewer, `--check` to
  validate/lint via the `ara` CLI).
- `/research-foresight <ara-dir> "<question>"` — to answer "what should I try
  next / why did this work / what if I change X", grounded in the artifact.
- `/submit-ara <dir>` — when an artifact is ready to publish to GitHub and list
  on the ARA Hub.

## Submission scope

When the researcher says “提交”, “重新提交”, or asks to push the current work,
inspect and publish all three repositories with native Git:

1. `/home/eric/Workspace/software_robot_library`
2. `/home/eric/Workspace/ws_trajectory_tracking/src`
3. `/home/eric/Workspace/MPCC_ARA`

Commit only validated, in-scope changes; do not create empty commits in clean
repositories. Push each current branch to its configured `origin`. Do not open
a pull request unless the researcher explicitly asks for one.

## Experiment record publication

After an experiment or numerical audit finishes and its results have been
validated, commit and push its canonical `/home/eric/Workspace/MPCC_ARA`
record automatically as part of the same task. Include the machine-readable
evidence and session trace. Do not wait for a separate submission request.

Do not publish a failed or still-unverified result as complete, and do not
bundle unrelated dirty-worktree changes. If the target branch, remote, or
in-scope file set is ambiguous, report that ambiguity before publishing.
