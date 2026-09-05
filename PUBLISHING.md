# Published on Zenodo

**Record:** <https://zenodo.org/records/22401515>
**Cite all versions:** [10.5281/zenodo.22401514](https://doi.org/10.5281/zenodo.22401514)
**This version (2.0.0):** [10.5281/zenodo.22401515](https://doi.org/10.5281/zenodo.22401515)

Published 5 September 2026 as an open-access preprint under CC BY 4.0, with the 75-page paper
and the 2.4 MB research package deposited together. Files are now immutable; metadata can still
be edited on Zenodo, and a new version can be added later without losing the concept DOI.

The checklist that produced it is kept below for the next release.

---

# Publishing checklist — Zenodo

Everything below is prepared. The steps that require signing in to an account are
yours to do; nothing in this repository logs in or uploads on your behalf.

## What is ready

| Item | Where | Status |
|---|---|---|
| Preprint PDF, 75 pages | `version2/output/pdf/Genetics_Complete.pdf` | ready |
| Markdown source | `version2/manuscript/genetics_unified.md` | ready |
| Full research package, 187 files | `version2/output/Genetics_Research_Package.zip` | ready, 2.4 MB |
| Zenodo metadata | `.zenodo.json` | ready |
| Citation metadata | `CITATION.cff` | ready |
| Licence, dual | `LICENSE` | CC BY 4.0 (writing) + MIT (code) |
| Availability, disclosure, scope | paper front matter | in the PDF |

The paper's front matter already carries the code and data availability statement,
the external-data provenance, a no-human-subjects statement, the AI-assistance
disclosure, the scope of claims and the licence.

## Route A — manual upload (simplest, works today)

1. Go to <https://zenodo.org> and sign in. You can sign in with GitHub or ORCID;
   no institutional affiliation is required.
2. **New upload**.
3. Drag in `Genetics_Complete.pdf`. Optionally also drag in
   `Genetics_Research_Package.zip` so the code, data and results are archived with
   the paper under the same DOI. It is 2.4 MB and contains exactly one PDF.
4. Fill the form using `.zenodo.json` as the source of truth:
   - **Upload type**: Publication → Preprint
   - **Title**: copy the `title` field
   - **Authors**: Mehran, Arshyia. Add your ORCID if you have one; if not, create
     one free at <https://orcid.org> first — it takes about two minutes and makes
     the record permanently attributable to you.
   - **Description**: copy the `description` field (it is already HTML)
   - **Licence**: Creative Commons Attribution 4.0 International
   - **Keywords**: copy the `keywords` list
   - **Related identifiers**: `https://github.com/Arshyi/GeneticInheritanceModelling`,
     relation "is supplement to"
5. **Publish**. Zenodo mints the DOI immediately.
6. Come back and add the DOI badge to `README.md` and the `doi:` field to
   `CITATION.cff`, then commit.

## Route B — GitHub release, auto-archived

This links the repository to Zenodo so every future release is archived
automatically. **The order matters**: enabling the repository must happen before
the release is created, or Zenodo will not see it.

1. Sign in to Zenodo with GitHub.
2. Go to <https://zenodo.org/account/settings/github/> and switch this repository
   **on**.
3. Only then, create the release:

```bash
git tag -a v2.0.0 -m "Unified preprint, 75 pages"
git push origin v2.0.0
```

4. On GitHub, publish a release from that tag. Zenodo picks it up within a minute
   or two and mints a DOI from `.zenodo.json`.
5. Add the DOI badge it gives you to `README.md`.

Route B archives the repository as a source snapshot. If you want the PDF itself
to be the primary record, do Route A as well — most people do both, and Zenodo
handles them as separate records.

## Before you press publish

- [ ] Read the front matter once more. It states that this is a preprint, that no
      predictive accuracy is established for any trait, and that Part VI is
      derivation only. Those statements need to stay true.
- [ ] Decide whether to include your ORCID. Recommended.
- [ ] Confirm you are happy with CC BY 4.0. It permits commercial reuse with
      attribution and is effectively irreversible once people rely on it.
- [ ] Nothing else of yours ships. The 54-page Version I write-up and the Extended
      Essay are not in the repository and not in the release zip; only their
      mathematics, their MATLAB and a SHA-256 digest are carried forward inside the
      paper itself.

## After publishing

- ~~Add the DOI to `README.md` and `CITATION.cff`.~~ Done.
- If you later pursue arXiv (q-bio.PE or cs.DS), you will need an endorsement from
  an existing arXiv author in that category. Having a Zenodo DOI and a public
  repository first makes that request considerably easier to make.
