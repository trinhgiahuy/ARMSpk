# Compiler comparison on A64FX

## Cloning the repo
```
git clone --recurse-submodules https://gitlab.com/domke/a64fxCvC
```

## License
BSD-3-Clause

## HowTo use this framework (follow these instructions):
### Installation of dependencies:
- check ./inst/_init.sh for missing dependencies

### Compile each benchmark:
- for llvm12 without polly do `cd inst; for x in `ls *sh`; do sed -i -e 's/flto=full/flto=thin/g' -e 's/-mllvm -polly -mllvm -polly-vectorizer=polly//' $x; done; cd -`
- execute `./inst/*.sh [fujitrad | fujiclang | llvm12 | gnu]` for all files and the desired compiler (except those starting with underscore)
- note: some benchmarks might need external and/or download additional files which are not part of this repo due to license contraints
- note: a newer version of RIKEN's micro kernels can be accessed here: https://github.com/RIKEN-RCCS/fs2020-tapp-kernels/tree/main/tapp-kernels ; please integrate it yourself if needed

### Running each benchmark:
- execute `./run/*/test.sh [fujitrad | fujiclang | llvm12 | gnu]` to search for best number ranks and threads per node
- modify BESTCONF variable in ./conf/*
- execute `./run/*/best.sh [fujitrad | fujiclang | llvm12 | gnu]` for final performance runs
- investigate logs in ./log/<hostname>/[b|t]estrun/<bench>.log

## Citing this work
Publications re-using (parts of) this repo and/or citing our work should refer to:
>J. Domke, "A64FX – Your Compiler You Must Decide!," in Proceedings of the 2021 IEEE International Conference on Cluster Computing (CLUSTER), EAHPC Workshop, (Portland, Oregon, USA), IEEE Computer Society, Sept. 2021.

Bibtex Entry:
```bibtex
@inproceedings{domke_a64fx_2021,
	address = {Portland, Oregon, USA},
	title = {{A64FX} -- {Your} {Compiler} {You} {Must} {Decide}!},
	booktitle = {2021 {IEEE} {International} {Conference} on {Cluster} {Computing} ({CLUSTER}), {EAHPC} {Workshop}},
	publisher = {IEEE Computer Society},
	author = {Domke, Jens},
	year = {2021},
}
```

## Development Team
Authored and maintained by Jens Domke ([@contact](http://domke.gitlab.io/#contact))
