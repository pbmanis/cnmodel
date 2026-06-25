## cnmodel — Claude Code context


## Project overview
Python package for running biophysically realistic models of auditory nerve/cochlear nucleus neurons.
Includes ion channel mechanisms, synaptic mechanisms, information about channel densities,
network connectivity and more. 

## Package layout
```
cnmodel/
  mechanisms/       # neuron mod files 
  data/   # text tables that define cells, channels, populations, connectivity and synapses
  decorator/            # tools to insert ion channels into specific parts of a cell with specific densities.
  guis/              # datatable GUI (PyQt) - dead code?
  populations/      # setup for specific cell types
  protocols/      # basic protocols to run for current clamp, voltage clamp, synapses, etc.
  synapses/         # warpper for different kinds of synapses and postsynaptic density mechanisms
  util/   # collected support functions
```

## Energy Use
- Do not explore large spaces in the solution landscape just because it is possible. Focus looking for solutions efficiently.
- 

## Coding conventions
- Never execute git commands.
- Never delete a file or directory.
- Never touch or modify files in the cnmodel/data directory
- Always ask before making any change to a source file, and make changes incrementally (ask for every change, even if they are related).
- When fixing a bug, comment out the original line(s) rather than deleting,
  and add a comment `# Claude fixed YYYY-MM-DD: <reason>` on the replacement.
- Use American, not British spellings.


## Testing / running
- No automated test suite covers the data pipeline end-to-end.
- "python test.py" tests specific key parts of the code and should be run after any change that might affect simulation results.
- `uv` is used for dependency management (`uv.lock`, `pyproject.toml`).
- Python 3.13 is required (pinned in `pyproject.toml`).

## Testing convention — REQUIRED after code changes
After editing any file in the following directories, run the relevant tests
before reporting a task complete:

  Directories:
  - cells
  - an_model
  - mechanisms
  - protocols
  - synapses
  - util (and its subdirectores)
  
  Commands:
    - Full suite (all of the above):
        `python test.py`

  Rules:
    - NEVER pass `--audit` to pytest; audit mode is for human use only
      at the terminal.
    - If any test fails, investigate and fix before declaring done.
  

## Common pitfalls
- Idioms from earlier versions of Neuron (before V9.0) may cause problems. Pay attention to the changes if tests fail.
- Comparing strings may fail because of case differences. Point these out.
- Early termination of a subroutine due to an error may not be properly propagated back to the top level, leading to subtle errors.
  
## Debugging
- Test for errors in logical conditions first (program flow, early termination, completeness, ambiguous conditions).
- Examine calculations for potential floating point errors.
- Rely on primary documentation for imported modules, not discussion boards.
- When looking for bugs/problems, and suggesting tests, explain why examining and testing specific scripts/sources is necessary.
- When potential bugs or fixes are found, the first action should always be to explain and ask before edits. 
- When debugging GUI-based interactions, always check the GUI for correctness first, before examining underlying computations.
- When debugging processes that file in a parallelization, it may be necessary to test a non-parallel version of the computation to pinpoint the actual error.
  -