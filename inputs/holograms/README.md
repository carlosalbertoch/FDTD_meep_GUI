# Input holograms

One folder per view, exactly as written by the acquisition system:

```
inputs/holograms/Phantom_Bottom_250/   i000.pgm … i249.pgm, Intensite_ref.pgm, config_manip.txt, recon.txt
inputs/holograms/Phantom_Top_250/      …
inputs/holograms/phantom_4pi/bottom/   32-hologram subset (used by inputs/config/phantom_4pi_subset.json)
inputs/holograms/phantom_4pi/top/      …
```

The folders are set in `views` of the configuration (`inputs/config/*.json`) or in the GUI
(*Input holograms*). Holograms are not uploaded to the repository (`.gitignore`).
