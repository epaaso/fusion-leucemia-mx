# Shared STAR Alignment for Arriba and STAR-Fusion

## Summary

Implement a new `--method arriba_starfusion` workflow that runs STAR once per
sample, streams BAM directly into Arriba, and retains
`Chimeric.out.junction` in the Nextflow work cache for concurrent STAR-Fusion
post-processing.

Validate a shared CTAT pseudo-masked GENCODE v50 reference first. If this
reduces clinically relevant Arriba sensitivity, especially `DUX4::IGH`, use
separate unmasked Arriba and pseudo-masked STAR-Fusion alignments.

## Reference Strategy

- Build and lock a GRCh38/GENCODE v50 CTAT reference bundle.
- Build a matching unmasked GRCh38/GENCODE v50 STAR index for the fallback
  Arriba workflow.
- Record genome assembly, GENCODE version, masking strategy, STAR version, and
  reference paths in logs and README.
- Never mix reference versions within a production cohort.

## Shared Workflow

Add `--method arriba_starfusion` while preserving existing `arriba`,
`starfusion`, and `fusioncatcher` modes.

```text
FASTQs
  |
  v
STAR_SHARED_ARRIBA
  |-- BAM stdout --> Arriba concurrently
  |
  `-- Chimeric.out.junction --> STARFUSION_CALL
```

### `STAR_SHARED_ARRIBA`

- Use the CTAT pseudo-masked STAR index.
- Preserve STAR-to-Arriba streaming through `--outStd BAM_Unsorted`.
- Produce both chimeric output forms:

```bash
--chimOutType Junctions WithinBAM HardClip
```

- Use the union of required sensitivity settings:

```bash
--outFilterMultimapNmax 50
--chimSegmentMin 10
--chimJunctionOverhangMin 10
--chimMultimapNmax 50
--chimOutJunctionFormat 1
--alignSJstitchMismatchNmax 5 -1 5 5
--peOverlapNbasesMin 10
```

- Emit:
  - Arriba fusion and discarded-fusion TSVs.
  - `Chimeric.out.junction` channel for STAR-Fusion.
  - `ReadsPerGene.out.tab`.
  - STAR logs.
- Keep junction and temporary alignment files only in the Nextflow work cache.

### `STARFUSION_CALL`

Run STAR-Fusion without realigning FASTQs:

```bash
STAR-Fusion \
  --genome_lib_dir <ctat_library> \
  --chimeric_junction Chimeric.out.junction \
  --output_dir starfusion_out
```

Emit normal STAR-Fusion predictions, abridged predictions, and logs.

## Efficient Scheduling

Configure the shared alignment and streamed Arriba process as the only
memory-heavy task:

```text
STAR_SHARED_ARRIBA:
  cpus: 16
  memory: 100.GB
  maxForks: 1
```

Configure post-alignment STAR-Fusion calls as lightweight concurrent tasks:

```text
STARFUSION_CALL:
  cpus: 1
  memory: 12.GB
  maxForks: 4
```

Expected behavior on the 20-CPU/125-GiB host:

- Only one STAR alignment runs at a time.
- Arriba runs concurrently with STAR through streaming.
- Up to four completed samples can undergo STAR-Fusion post-processing
  concurrently.
- STAR-Fusion calls may overlap with the next STAR/Arriba alignment.
- Do not allocate extra CPUs to `STARFUSION_CALL`; its standard
  junction-processing stages are single-threaded.
- Benchmark 16 versus 19 STAR alignment threads. Retain 16 only if total cohort
  throughput improves without excessive alignment slowdown.

## Validation and Fallback

Validate using:

- `LAL1276_S9`: known `DUX4::IGH` candidate.
- `103-MO-15-563486477`: strong `TCF3::PBX1` control.

Compare shared-workflow outputs against current caller-specific outputs:

- Clinically relevant fusion presence.
- Gene partners and breakpoint coordinates.
- Arriba confidence and supporting-read counts.
- STAR-Fusion junction and spanning-fragment counts.
- Total calls and unexpected high-confidence call changes.
- Runtime, peak memory, CPU utilization, and disk usage.

Shared validation fails if:

- Expected `DUX4::IGH` or `TCF3::PBX1` evidence disappears.
- A clinically relevant fusion is lost.
- Support decreases materially without a defensible explanation.
- Reference masking causes unacceptable Arriba behavior.

If validation fails, production uses:

```text
FASTQs
|-- STAR unmasked v50 | Arriba streaming
`-- STAR CTAT pseudo-masked v50 -> STAR-Fusion
```

Fallback scheduling:

- Run only one memory-heavy alignment at a time.
- Preserve STAR-to-Arriba streaming.
- Run lightweight STAR-Fusion post-processing concurrently when resources
  permit.

## Implementation Changes

- Refactor `nextflow_pipeline/main.nf`:
  - Add `arriba_starfusion` method validation and routing.
  - Add `STAR_SHARED_ARRIBA`.
  - Add `STARFUSION_CALL`.
  - Route caller outputs independently into `SUMMARIZE`.
  - Generate separate Arriba and STAR-Fusion cohort summaries.
- Extend `nextflow_pipeline/nextflow.config`:
  - Shared reference and resource parameters.
  - `shared_star_threads = 16`.
  - `shared_star_ram = 100.GB`.
  - `starfusion_call_threads = 1`.
  - `starfusion_call_ram = 12.GB`.
  - `starfusion_call_max_forks = 4`.
- Update `README.md`:
  - Document combined mode, reference locking, masking behavior, scheduling,
    validation, and fallback workflow.
- Preserve existing caller-specific modes and outputs unchanged.
- Do not alter or invalidate currently running cohort workflows.

## Tests

- Confirm STAR generates both streamed BAM and `Chimeric.out.junction`.
- Confirm Arriba consumes BAM from stdin while STAR is running.
- Confirm STAR-Fusion accepts the junction file and does not launch STAR.
- Confirm Nextflow starts STAR-Fusion calls as soon as each shared alignment
  finishes.
- Confirm multiple STAR-Fusion calls run concurrently while only one STAR
  alignment runs.
- Confirm `-resume` reuses shared alignment and caller outputs.
- Confirm work-cache-only intermediate retention.
- Confirm separate caller summaries contain all processed samples.
- Execute the two-sample validation before enabling combined mode for the full
  cohort.

## Assumptions

- Shared CTAT pseudo-masked alignment is experimental until validation passes.
- Dual alignment remains the accepted fallback.
- GENCODE v50 becomes the locked reference version only after validation.
- Existing v41/v44 outputs remain historical and are not silently mixed with
  new v50 cohort results.
