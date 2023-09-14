06/14/22

1) Tritrig_beam and wab_beam recon samples already exist: https://confluence.slac.stanford.edu/display/hpsg/pass0+for+3pt7+GeV
    Readout steering-file: /org/hps/steering/readout/PhysicsRun2019TrigMultiSingles.lcsim
    Recon steering-file: UNKNOWN-WAITING FOR TONGTONG

2) Need to generate Rad and mix with beam 
    Make sure target = 8um and beam current is adjusted in hps-mc run_parameters
    change nelec to 625  DONE

3) Need to generate SIMPs
    Make sure target = 8 um and beam current is adjusted to 50 nA in hps-mc DONE

 4) Use recon steering PhysicsRun2019MCRecon.lcsim with beamZ at 0.0

