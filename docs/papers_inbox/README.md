# Paper inbox

Drop papers here under whatever name they downloaded with. Everything except this README is
git-ignored, so the directory itself survives in the repo while its contents never do.

**Do not file anything here by hand.** Claude identifies each file from its own embedded
metadata and first page, renames it to the convention in [../papers/README.md](../papers/README.md),
moves it to `docs/papers/`, and adds a manifest row. Publisher filenames are not trustworthy —
Elsevier `mmc` numbers have twice turned out to be main articles rather than supplements — so
identification is always by content, never by the name it arrived with.

`scripts/check_papers.py` reports when files are sitting here waiting to be filed.

An empty inbox is the normal state.
