spinham/
├─── ...
├─── cell/ (attr: units)
│    ├──> a1
│    ├──> a2
│    └──> a3
├─── atoms/
│    ├─── 1/
│    │    ├──> name
│    │    ├──> position  # Always relative to the cell
│    │    ├──> spin
│    │    ├──> g-factor
│    │    ├──> charge
│    │    ├──> orbital-moment
│    │    └──> total-moment
│    ├─── ...
│    └─── n/
│         ├──> name
│         ├──> position  # Always relative to the cell
│         ├──> spin
│         ├──> g-factor
│         ├──> charge
│         ├──> orbital-moment
│         └──> total-moment
├─── notation/
│    ├──> double-counting
│    ├──> spin-normalized
│    ├──> exchange-factor
│    └──> on-site-factor
├─── exchange/ (attr: units)
│    ├─── 1/
│    │    ├──> atom-1
│    │    ├──> atom-2
│    │    ├──> ijk
│    │    └──> matrix
│    ├─── ...
│    └─── n/
│         ├──> atom-1
│         ├──> atom-2
│         ├──> ijk
│         └──> matrix
├─── on-site/ (attr: units)
│    ├─── 1/
│    │    ├──> atom
│    │    └──> matrix
│    ├─── ...
│    └─── n/
│          ├──> atom
│          └──> matrix
└─── ...
