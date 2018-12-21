# SN_MAF
SN_MAF is a framework for estimating MAF metrics for supernovae (SN). It can also be used to simulate SN light curves. It proposes a set of python script that can be run along with yaml input files where specific parameters may be chosen. 


- **Instruction for installation**
  - git clone https://github.com/lsstdesc/SN_MAF.git
  - cd SN_MAF
  - git checkout dev

- **Running SN-MAF metrics**

  - Environnement setup
    - setup script: SN_MAF/setups/setup_metric.sh
    - This script requires up to two arguments and can be run:
      - @NERSC: source SN_MAF/setups/setup_metric.sh NERSC
      - @CCIN2P3: source SN_MAF/setups/setup_metric.sh CCIN2P3
    - if you wish to run elsewhere then you need to provide the full path to the setup script corresponding to a release including the lsst_sims package, ie source SN_MAF/setups/setup_metric.sh MYENV full_path_to_setup_script_stack.


  - Running the Cadence metric
    - python SN_MAF/run_scripts/run_cadence_metric.py SN_MAF/input/param_cadence_metric.yaml
    - A description of the input yaml file is given [here](doc/yaml_cadence.md)
    - you may have to change the 'filename' parameter to the OpSim db name you would like to use as input.
    - output : a set of plots: 
      - Mean cadence vs mean m5 (5-sigma depth) <img src="doc/cadence_m5_r.png" height="24">
      - Histogram of redshift limits <img src="doc/zlim_r.png" height="24">

  - Running the Signal-to-Noise Ratio (SNR) metric
    -  python SN_MAF/run_scripts/run_snr_metric.py SN_MAF/input/param_snr_metric.yaml
    - A description of the input yaml file is given [here](doc/yaml_snr.md)
    - you may have to change the 'filename' parameter to the OpSim db name you would like to use as input.
    - output : a set of plots:
      - SNR vs Time (per band and per season) <img src="doc/snr_z_season_1.png" height="24">