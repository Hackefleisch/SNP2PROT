# Domain policy — what `dbd_seq` is, and what is excluded

> **`dbd_seq` is the Pfam envelope padded by 10 residues on each side. It is not the bare
> Pfam domain, and it is not the sequence that was physically on the array.**
> The padding is deliberate. Read §2 before using the column for anything.

Machine-readable settings: the `domain:` block in
[../configs/thresholds.yaml](../configs/thresholds.yaml). Implementation:
[../src/snp2prot/domains.py](../src/snp2prot/domains.py).

## 1. The three conditions

A construct is admitted only if the stored subunit satisfies all three. They come from what
the dataset is *for*: predicting a structure, locating the centre of the DNA-contacting
residues, and embedding the residues around that centre.

1. **Sole responsibility.** The stored subunit alone produced the measured interaction.
2. **One continuous region.** A single contiguous stretch, not fragments scattered through
   the protein — otherwise "the centre of the contact residues" is a point in space whose
   surroundings are an artefact of how the pieces were arranged.
3. **The variation lies inside it.** Every change that alters binding behaviour falls within
   the stored region. A mutation outside it makes a variant sequence-identical to its own
   wild type while carrying a different label — invisible to a sequence model and to a
   structure model alike.

Constructs that fail are **dropped and reported**, never silently repaired.

## 2. Why the padding exists

Pfam's models begin at the structural core. For a homeodomain that is the helix-turn-helix,
which clips the **N-terminal arm** — and the N-terminal arm reaches into the DNA minor groove
and makes genuine base contacts.

This is not a theoretical worry. Of BAR15A's 91 point mutations in single-domain constructs,
**four fall outside the bare Pfam envelope**, three of them N-terminal by 2–9 residues:

| variant | mutation at | Pfam envelope | offset | trimmed to the bare envelope |
|---|---|---|---|---|
| `KLF11_R402Q` | 13 | zf-C2H2 [15–39] | 2 before | identical to `KLF11_REF` |
| `VSX1_G160D` | 10 | Homeodomain [15–71] | 5 before | identical to `VSX1_REF` |
| `SNAI2_D119E` | 6 | zf-C2H2 [15–37] | 9 before | identical to `SNAI2_REF` |
| `ZNF655_E327G` | 137 | zf-C2H2 [113–135] | 2 after | inside the array span |

`padding_aa: 10` covers all of them with margin. The value is set from the biology of the
N-terminal arm, not fitted to these four cases; `domains.audit_variant_positions` re-checks
condition 3 on every variant at parse time, so a future source that violates it is rejected
rather than quietly flattened.

**Consequence for downstream use:** `dbd_seq` contains ~10 residues of flanking sequence on
each side of the Pfam domain. If you need the bare domain, the protein table carries the
unpadded envelope offsets.

## 3. What is excluded, and why

| rejection | condition | what it is |
|---|---|---|
| `mixed_families` | 1, 2 | Pfam hits from two different families in one construct — PAX + homeodomain, POU + homeodomain. Two domains bind; nothing in the data says which produced a given score, and their relative geometry is undetermined without DNA. |
| `repeat_array` | 2 | Repeated hits of one family — a C2H2 zinc-finger array is 2–6 separate ~23-residue folds on flexible linkers, each requiring its own Zn²⁺. Not one continuous unit. |
| `no_domain` | — | No Pfam hit above the family's gathering threshold. We cannot say what is doing the binding. |
| `unresolved_residue` | — | The padded window contains a character outside the 20 standard residues. `X`, `B`, `Z` and `U` record the depositor's uncertainty, not the protein: a sequence model would embed the ambiguity as a token and a structure predictor cannot place the residue at all. One domain in the corpus was affected, `MAR17A:Esrrb` with two `X` in an 89 aa zf-C4 domain, excluded 2026-08-18 (`docs/DECISIONS.md`, `D2`). |
| mutation outside the padded domain | 3 | The variant would be sequence-identical to its wild type. |
| protein complex | 1 | Two or three *different chains* forming one binding unit — `Myc_Max`, `Kay_Jra` (Fos/Jun), `Da_Twi`, the C. elegans `HLH-2_*` heterodimers, the `MAML1-CSL-GST-NOTCH*` ternary complexes. No single `dbd_seq` is responsible for the measurement. This is the same principle as the mixed-family rejection, one level up: there the two domains sat in one chain, here they sit in two. UniPROBE publishes no sequence for any of them, so they are excluded either way, but the reason is now recorded rather than incidental. |

## 3b. One domain described by several Pfam models

Pfam sometimes models a single domain with more than one family. Applied naively, the
`mixed_families` rule then rejects a perfectly good protein, because two family names appear
in one construct.

This is not hypothetical. An audit of every construct found **12 rejections where the two
families covered the same region**, eight of them bZIP:

| construct | hits | overlap |
|---|---|---|
| `MAR17A:Atf3` | `bZIP_1[18-74]` + `bZIP_2[18-70]` | 100% |
| `SCI09:Jundm2` | `bZIP_1[16-74]` + `bZIP_2[17-70]` | 100% |
| `SHO18A:Jra` | `bZIP_1[210-273]` + `bZIP_2[216-264]` | 100% |
| `SCI09:Zfp691` | `zf-H2C2_2[30-55]` + `zf-C2H2[44-66]` | 52% |
| `EMBO10:Gm4881` | `Ets[16-96]` + `HSF_DNA-bind[23-88]` | 100% |

bZIP is one of the three families the brief names explicitly, and it was being discarded
almost entirely, with rejection messages that read as if the policy were working.

**Rule.** Before the policy is applied, hits overlapping more than 50% of the shorter one are
collapsed to the highest-scoring, since they describe one domain. Deliberately narrow:

- separated hits still reject — PAX + homeodomain remain two spatial units (condition 1);
- partially overlapping hits still reject — two domains touching is not one domain twice;
- zinc-finger arrays remain arrays — fingers sit side by side, so `repeat_array` still fires.

Recovering bZIP took it from 6 to 13 domains.

## 3c. Family naming

Collapsing leaves whichever model scored higher, so `dbd_family` would otherwise depend on an
implementation detail: bZIP appeared as `bZIP_1` for some proteins and `bZIP_2` for others,
splitting one family in two. Since `dbd_family` drives leave-one-family-out, a family split
this way could never be held out as a group.

Synonymous model names are therefore mapped to one label (`domains.FAMILY_ALIASES`):
`bZIP_1`/`bZIP_2` → `bZIP`, `zf-H2C2_2` → `zf-C2H2`, `Homeobox_KN` → `Homeodomain`,
`SOXp` → `HMG_box`.

**`TF_AP-2` (mammalian TFAP2) and `AP2` (plant AP2/ERF) are deliberately NOT merged** — similar
names, different domains. A test pins this so the two are not tidied together later.

## 4. Consequences for the phase plan

**Phase 3 (Persikov B1H + Najafabadi C2H2) is dropped.** Both are C2H2 zinc-finger arrays, so
every entry fails condition 2 — roughly 8,000 domains that cannot be admitted. Persikov's
design makes it sharper still: the variable finger sits inside a fixed three-finger context,
so the subunit that *varies* and the subunit that *binds* are deliberately different, which
fails conditions 1 and 3 simultaneously. Its selection-absence negatives were already a
concern (open item #1); this settles it.

**Phases 4 and 6 need per-construct screening, not wholesale acceptance or rejection.**
SNP-SELEX (Phase 4) varies the DNA rather than the protein, so condition 3 is vacuous on the
protein axis — but many of its 270 TFs are C2H2 and will fail condition 2. The Tier 4 test
sets (Phase 6) are bHLH proteins (Pho4, Cbf1, MAX), which bind as dimers; MAX's substitutions
are described as being "in and around" the DBD, so condition 3 has to be checked per variant.
Each is screened by the same policy when its phase is built.

## 5. Cluster membership and indels

Members of a `wt_id` cluster need **not** be the same length. Distance comes from a pairwise
alignment (`snp2prot.align`), so a variant may carry insertions or deletions.

Terminal gaps are free, and that is not a convenience — it is required for correctness here.
`dbd_seq` is the padded envelope clipped where the assayed construct ends, so the amount of
padding that survives depends on the construct, not on the protein:

| ROG18A construct | bare Pfam domain | padded `dbd_seq` | N-terminal padding kept |
|---|---:|---:|---:|
| `FoxJ3` | 83 aa | 103 aa | 10 |
| `FoxJ3_N3_6aa` | 83 aa | 95 aa | 2 |

Those two differ by 8 residues at the terminus and by **zero** residues of biology. Charging
for the terminal gap would report 8 phantom indels between a protein and its own chimera.
The validator now *warns* rather than errors on a length-heterogeneous cluster, and the
warning says to check whether the difference is a real indel or clipped padding.

## 6. Changing the policy

Every knob is in the `domain:` block of `configs/thresholds.yaml`. Setting
`allow_repeat_arrays: true` readmits the C2H2 arrays as a single spanned region;
`allow_mixed_families: true` readmits PAX and POU. Both change what the dataset *is*, so
they invalidate every table, report and provenance row — rebuild from `data/raw/`, which is
untouched by any of this.
