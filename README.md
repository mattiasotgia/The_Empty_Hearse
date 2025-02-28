# Readme intro

This is (starting from today, feb. 26th, 2025) the main folder where I will try to keep all the tidy documents 

# Selection cuts (MA selection)

At document top there are the standard definitions (these are in `namespace maria_selection`)

 - `isInFV` checks if (x, y, z) point is in fiducial volume (with some sanity checks on their values and accounting for dangling wire)
 - `isInContained` 5 cm containment for Y in both cryo $\to$ set `dist = 5` (default value set)
 - `isInActive` check if is in the active TPC volume
 - `isInDetector` equivalent to `isInContained` with `dist=5`.

These are used in the automatic selections

 - `all_contained` take a `const caf::Proxy<caf::SRSlice> slice` and run through `std::vector<caf::SRPFP> slice.reco.pfp`. 