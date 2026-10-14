# sst-milne — notes for agents

Milne-Eddington inversions (pyMilne) of SST/CRISP Fe I 6302 SSTRED cubes,
with quicklook plots/movies and Slurm scripts.

- `README.md`: installation (Mac, Linux, Slurm) and usage.
- `HANDOVER.md`: where data/results live, run history, problems already
  solved (read before debugging), open decisions.

Working conventions with the owner (Jonas Zbinden):
- He runs jobs on the sml server and the cluster himself; give him
  complete, copy-paste commands (one per code block) and ask for the log
  lines (`grep '^scan'`) or figures back.
- Never commit pyMilne sources here or modify the upstream pyMilne
  checkout; machine-specific build changes go into `pymilne_build/`.
- Test every change locally on a test cube (`tools/make_test_cube.py`)
  with a short `--scans 0-1 --crop ...` run before pushing.
- Results (FITS, movies) never go into git (`.gitignore`).
- Commit messages explain why; end them with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Ask before deleting results or anything on the machines.
