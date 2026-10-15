# examples/

* `user_paths_handoff_layout.json` — a filled-in `config/user_paths.json` for
  the layout of the extracted full handoff archive (repository and `data/`
  side by side). Copy it to `config/user_paths.json` and adjust if your data
  lives elsewhere. Derived from the handoff's generated
  `PORTABLE_USER_PATHS_SUGGESTION.json`.

Setup is always the same: edit **one** paths file (`config/user_paths.json`),
keep `config/active_profile.json` at `"profile": "user"`, then run
`notebooks/00_Explore_Data.ipynb` onward in order.
