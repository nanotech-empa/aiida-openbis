```mermaid
erDiagram
Container {

}
WaferSubstrate {
    string perm_id  
    string name  
    string batch  
    string sample_plate  
    boolean hazardous  
    string hazardous_specification  
    boolean fridge  
    boolean no_light  
    boolean dry  
    boolean no_oxygen  
    boolean other_storage_condition  
    string other_storage_condition_specification  
    date receive_date  
    string comments  
}
Chemist {
    string perm_id  
    string name  
}
Supplier {
    string perm_id  
    string name  
    string email  
    string work_phone  
    string address  
}
Person {
    string perm_id  
    string name  
    string short_name  
    boolean work_status  
    string email  
    string work_phone  
    string mobile_phone  
}
Storage {
    string perm_id  
    string name  
    string comments  
}
Room {
    string perm_id  
    string name  
    string comments  
}
Institution {
    string perm_id  
    string name  
    string description  
    string address  
}
QuantityValue {
    double has_value  
    string has_unit  
}
SubstratesDetails {
    string material  
    string dopants  
    string coverage  
}
UHVComponent {
    string perm_id  
    string name  
    string description  
    UHVComponentTypeEnum uhv_component_type  
    string empa_id  
    string model  
    string serial_number  
    date receive_date  
    string comments  
}
Manufacturer {
    string perm_id  
    string name  
    string email  
    string work_phone  
}
Transfer {
    string perm_id  
    string name  
    string description  
}
Sample {
    string perm_id  
    string name  
}
Layered2DMaterial {
    string perm_id  
    string name  
    integer empa_number  
    string batch  
    string sample_plate  
    boolean hazardous  
    string hazardous_specification  
    boolean fridge  
    boolean no_light  
    boolean dry  
    boolean no_oxygen  
    boolean other_storage_condition  
    string other_storage_condition_specification  
    date receive_date  
    string supplier_own_name  
    string comments  
}
Layers2DDetails {
    string material  
    integer number_layers  
    string dopants  
    string coverage  
}
Crystal {
    string perm_id  
    string name  
    string material  
    string face  
    string sample_plate  
    boolean hazardous  
    string hazardous_specification  
    boolean fridge  
    boolean no_light  
    boolean dry  
    boolean no_oxygen  
    boolean other_storage_condition  
    string other_storage_condition_specification  
    string specifications  
    date receive_date  
    string reference_number  
    string comments  
}
TipSensor {
    string perm_id  
    string name  
    TipSensorTypeEnum tip_sensor_type  
    string wire_material  
    string holder_id  
    double q_factor  
    string reference_number  
    string device_parameters  
    date receive_date  
    string comments  
}
TipPreparation {
    string perm_id  
    string name  
    string description  
    TipPreparationTypeEnum tip_preparation_type  
}
Author {
    string perm_id  
    string name  
}
TightBinding {
    string perm_id  
    string name  
    string wfms_url  
    string wfms_uuid  
    stringList hoppings  
    string lattice_type  
    stringList hamiltonian  
    string comments  
}
KPointsConditions {

}
AtomsPositions {
    string atom_symbol  
}
OneDMeasurement {
    string perm_id  
    string name  
    custom_datetime start_time  
    ExperimentTypeEnum experiment_type  
    FilterTypeEnum filter_type  
    integer filter_order  
    integer num_pixel  
    boolean backward_sweep  
    integer num_sweeps  
    stringList channel_names  
    boolean record_final_z  
    boolean bias_spectroscopy_reset_bias  
    boolean bias_spectroscopy_z_controller_hold  
    boolean z_spectroscopy_reset_z  
    boolean feedback_active  
    FeedbackTypeEnum feedback_type  
    double z_controller_p_gain  
    double z_controller_i_gain  
    boolean lock_in_status  
    string lock_in_modulated_signal  
    integerList lock_in_hp_filter_orders  
    LockInDemodulatedSignalEnum lock_in_demodulated_signal  
    integerList lock_in_harmonics  
    integerList lock_in_lp_filter_orders  
    boolean lock_in_sync_filter  
    string wfms_uuid  
    string comments  
}
TwoDMeasurement {
    string perm_id  
    string name  
    custom_datetime start_time  
    boolean feedback_active  
    FeedbackTypeEnum feedback_type  
    double z_controller_p_gain  
    double z_controller_i_gain  
    string comments  
    PModelEnum p_model  
    string wfms_uuid  
}
Instrument {
    string perm_id  
    string name  
    string serial_number  
    string empa_id  
    string documentation_website  
    date receive_date  
    string comments  
}
DFT {
    string perm_id  
    string name  
    string wfms_uuid  
    string wfms_url  
    double fermi_energy  
    XCFunctionalEnum xc_functional  
    boolean spin_polarised_calculation  
    doubleList initial_spin_guesses  
    boolean spin_multiplicity  
    double total_spin_squared  
    double net_charge  
    boolean force_multiplicity  
    string comments  
}
OutputPopulationAnalysis {
    double charge  
    double magnetic_momentum  
}
Deposition {
    string perm_id  
    string name  
    string comments  
}
Molecule {
    string perm_id  
    string name  
    string iupac_name  
    string sum_formula  
    string smiles  
    string cas_number  
    integer empa_number  
    string batch  
    string vial  
    boolean hazardous  
    string hazardous_specification  
    boolean fridge  
    boolean no_light  
    boolean dry  
    boolean no_oxygen  
    boolean other_storage_condition  
    string other_storage_condition_specification  
    string chemist_molecule_name  
    date receive_date  
    string comments  
}
EvaporationTemperature {
    string device  
    custom_datetime datetime  
}
Annealing {
    string perm_id  
    string name  
    string comments  
}
Sputtering {
    string perm_id  
    string name  
    string comments  
}
EvaporatorSlot {
    integer evaporator_number  
    string details  
}
OscillationControlSettings {
    boolean differential_input  
    boolean one_ten  
}
PLLSetupSettings {
    double q_factor  
}
OutputSettings {
    boolean output_active  
    boolean output_add  
}
AmplitudeSettings {
    boolean controller_on  
}
PhaseSettings {
    boolean controller_on  
}
DemodulationSettings {
    doubleList inputs  
    doubleList filters_orders  
}
ScanSettings {
    doubleList field  
    DirectionSpeedEnum direction_speed  
}
ScanSpeed {

}
ScanPixels {

}
ScanRange {

}
ScanTime {

}
ScanOffset {

}
PiezoConfigurationSettings {
    ActiveCalibEnum active_calib  
    boolean drift_correction_status  
}
DriftSettings {

}
SecondOrderCorrection {

}
CurvatureRadius {

}
ScanSlope {

}
HVGainSettings {

}
SensitivitySettings {

}
SRD {
    string perm_id  
    string name  
    string comments  
}
Software {
    string perm_id  
    string name  
    string repository_url  
}
Result {
    string perm_id  
    string name  
    string description  
    string comments  
}
GeometryOptimisation {
    string perm_id  
    string name  
    string wfms_uuid  
    string wfms_url  
    doubleList geometry_constraints  
    string comments  
}
EmpiricalModelling {
    string perm_id  
    string name  
    string wfms_uuid  
    string wfms_url  
}
AtomisticModel {
    string perm_id  
    string name  
    string comments  
    boolean optimised  
    string wfms_uuid  
}
PeriodicBoundaryConditions {

}
ForceConvergenceThreshold {
    double charge  
    double magnetic_momentum  
}
Publication {
    string perm_id  
    string name  
    string abstract  
    string doi  
    integer year  
    string url  
    string dataset_url  
    string comments  
}
Draft {
    string perm_id  
    string name  
    DraftTypeEnum draft_type  
    string comments  
}
Grant {
    string perm_id  
    string name  
    string funder  
    date start_date  
    date end_date  
    string project_id  
    string acknowledgement_sentence  
    string comments  
}
Protocol {
    string perm_id  
    string name  
    string description  
}
PDOS {
    string perm_id  
    string name  
    string wfms_uuid  
    string wfms_url  
    doubleList atomic_selections  
    doubleList orbitals_selections  
    PlotContributionsEnum plot_contributions  
    PlotGroupByEnum plot_group_by  
    string comments  
}
OpticalComponent {
    string perm_id  
    string name  
}
Notes {
    string perm_id  
    string name  
    string description  
}
MinimumEnergyPotential {
    string perm_id  
    string name  
    string wfms_uuid  
    string wfms_url  
    MethodTypeEnum method_type  
    integer number_geometries  
    doubleList geometry_constraints  
    doubleList geometry_constraints_increments  
    string comments  
}
Maintenance {
    string perm_id  
    string name  
    string description  
    string comments  
}
HydrogenCracker {
    string perm_id  
    string name  
    string sum_formula  
    string comments  
}
FillCryostat {
    string perm_id  
    string name  
    string substance  
    DewarEnum dewar  
    string comments  
}
Errors {
    string perm_id  
    string name  
    string description  
}
Dosing {
    string perm_id  
    string name  
    string sum_formula  
    string comments  
}
CanGasBottle {
    string perm_id  
    string name  
    string description  
    string state  
    string comments  
}
BandStructure {
    string perm_id  
    string name  
    string wfms_uuid  
    string wfms_url  
    string comments  
}
BSEnergy {
    double charge  
    double magnetic_momentum  
}

Container ||--}o Annealing : "annealings"
Container ||--}o AtomisticModel : "atomisticmodels"
Container ||--}o Author : "authors"
Container ||--}o BandStructure : "band_structures"
Container ||--}o CanGasBottle : "can_gas_bottles"
Container ||--}o Chemist : "chemists"
Container ||--}o Crystal : "crystals"
Container ||--}o Deposition : "depositions"
Container ||--}o DFT : "dfts"
Container ||--}o Dosing : "dosings"
Container ||--}o Draft : "drafts"
Container ||--}o EmpiricalModelling : "empirical_modellings"
Container ||--}o Errors : "errors"
Container ||--}o FillCryostat : "fillcryostats"
Container ||--}o GeometryOptimisation : "geoopts"
Container ||--}o Grant : "grants"
Container ||--}o HydrogenCracker : "hydrogen_crackers"
Container ||--}o Institution : "institutions"
Container ||--}o Instrument : "instruments"
Container ||--}o Layered2DMaterial : "layered2dmaterials"
Container ||--}o Maintenance : "maintenances"
Container ||--}o Manufacturer : "manufacturers"
Container ||--}o MinimumEnergyPotential : "meps"
Container ||--}o Molecule : "molecules"
Container ||--}o Notes : "notes"
Container ||--}o OpticalComponent : "optics"
Container ||--}o PDOS : "pdoss"
Container ||--}o Person : "persons"
Container ||--}o Protocol : "protocols"
Container ||--}o Publication : "publications"
Container ||--}o Result : "results"
Container ||--}o Room : "rooms"
Container ||--}o Sample : "samples"
Container ||--}o Software : "softwares"
Container ||--}o Sputtering : "sputterings"
Container ||--}o SRD : "srds"
Container ||--}o TwoDMeasurement : "two_d_measurements"
Container ||--}o Storage : "storages"
Container ||--}o OneDMeasurement : "one_d_measurements"
Container ||--}o Supplier : "suppliers"
Container ||--}o TightBinding : "tight_bindings"
Container ||--}o TipPreparation : "tip_preparations"
Container ||--}o TipSensor : "tip_sensors"
Container ||--}o Transfer : "transfers"
Container ||--}o UHVComponent : "uhv_components"
Container ||--}o WaferSubstrate : "wafer_substrates"
WaferSubstrate ||--}o SubstratesDetails : "substrates"
WaferSubstrate ||--|o QuantityValue : "diameter"
WaferSubstrate ||--|o QuantityValue : "height"
WaferSubstrate ||--|o QuantityValue : "thickness"
WaferSubstrate ||--|o Storage : "storage"
WaferSubstrate ||--|o Chemist : "chemist"
Chemist ||--|o Person : "person"
Chemist ||--|o Supplier : "supplier"
Storage ||--|o Room : "room"
Room ||--|o Institution : "institution"
SubstratesDetails ||--|o QuantityValue : "thickness"
UHVComponent ||--|o Manufacturer : "manufacturer"
Transfer ||--|o Room : "location_before"
Transfer ||--|o Room : "location_after"
Transfer ||--|o Sample : "sample"
Sample ||--|o Crystal : "crystal"
Sample ||--|o Layered2DMaterial : "layered_2d_material"
Sample ||--|o WaferSubstrate : "wafer_substrate"
Layered2DMaterial ||--}o Layers2DDetails : "layers_2d"
Layered2DMaterial ||--}o SubstratesDetails : "substrates"
Layered2DMaterial ||--|o QuantityValue : "width"
Layered2DMaterial ||--|o QuantityValue : "length"
Layered2DMaterial ||--|o QuantityValue : "thickness"
Layered2DMaterial ||--|o Storage : "storage"
Layered2DMaterial ||--|o Chemist : "chemist"
Crystal ||--|o QuantityValue : "diameter"
Crystal ||--|o QuantityValue : "height"
Crystal ||--|o Storage : "storage"
Crystal ||--|o Supplier : "supplier"
TipSensor ||--|o QuantityValue : "thickness"
TipSensor ||--|o QuantityValue : "length"
TipSensor ||--|o QuantityValue : "resonance_frequency"
TipSensor ||--|o Manufacturer : "manufacturer"
TipPreparation ||--|o TipSensor : "tip_sensor"
TipPreparation ||--}o Author : "authors"
Author ||--|o Person : "person"
Author ||--|o Institution : "institution"
TightBinding ||--}o AtomsPositions : "atoms_positions"
TightBinding ||--}o QuantityValue : "cell_vectors"
TightBinding ||--|o KPointsConditions : "k_points_conditions"
KPointsConditions ||--|o QuantityValue : "x_coordinate"
KPointsConditions ||--|o QuantityValue : "y_coordinate"
KPointsConditions ||--|o QuantityValue : "z_coordinate"
AtomsPositions ||--|o QuantityValue : "x_coordinate"
AtomsPositions ||--|o QuantityValue : "y_coordinate"
AtomsPositions ||--|o QuantityValue : "z_coordinate"
OneDMeasurement ||--|o QuantityValue : "duration"
OneDMeasurement ||--|o QuantityValue : "acquisition_coordinate_x"
OneDMeasurement ||--|o QuantityValue : "acquisition_coordinate_y"
OneDMeasurement ||--|o QuantityValue : "acquisition_coordinate_z"
OneDMeasurement ||--|o QuantityValue : "final_z"
OneDMeasurement ||--|o QuantityValue : "filter_cutoff"
OneDMeasurement ||--|o QuantityValue : "bias_setpoint"
OneDMeasurement ||--|o QuantityValue : "bias_calibration_factor"
OneDMeasurement ||--|o QuantityValue : "bias_calibration_offset"
OneDMeasurement ||--|o QuantityValue : "z_avg_time"
OneDMeasurement ||--|o QuantityValue : "first_settling_time"
OneDMeasurement ||--|o QuantityValue : "settling_time"
OneDMeasurement ||--|o QuantityValue : "integration_time"
OneDMeasurement ||--|o QuantityValue : "end_settling_time"
OneDMeasurement ||--|o QuantityValue : "z_control_time"
OneDMeasurement ||--|o QuantityValue : "max_slew_rate"
OneDMeasurement ||--|o QuantityValue : "bias_spectroscopy_sweep_start"
OneDMeasurement ||--|o QuantityValue : "bias_spectroscopy_sweep_end"
OneDMeasurement ||--|o QuantityValue : "bias_spectroscopy_z_offset"
OneDMeasurement ||--|o QuantityValue : "z_spectroscopy_initial_z_offset"
OneDMeasurement ||--|o QuantityValue : "z_spectroscopy_sweep_distance"
OneDMeasurement ||--|o QuantityValue : "z_spectroscopy_time_between_forward_backward"
OneDMeasurement ||--|o QuantityValue : "current_setpoint"
OneDMeasurement ||--|o QuantityValue : "current_calibration_factor"
OneDMeasurement ||--|o QuantityValue : "current_calibration_offset"
OneDMeasurement ||--|o QuantityValue : "current_gain"
OneDMeasurement ||--|o QuantityValue : "z_position"
OneDMeasurement ||--|o QuantityValue : "z_controller_setpoint"
OneDMeasurement ||--|o QuantityValue : "z_controller_time_constant"
OneDMeasurement ||--|o QuantityValue : "z_controller_tip_lift"
OneDMeasurement ||--|o QuantityValue : "z_controller_switch_off_delay"
OneDMeasurement ||--|o QuantityValue : "lock_in_frequency"
OneDMeasurement ||--|o QuantityValue : "lock_in_amplitude"
OneDMeasurement ||--}o QuantityValue : "lock_in_hp_filter_cutoffs"
OneDMeasurement ||--}o QuantityValue : "lock_in_reference_phases"
OneDMeasurement ||--}o QuantityValue : "lock_in_lp_filter_cutoffs"
OneDMeasurement ||--|o QuantityValue : "sample_temperature"
OneDMeasurement ||--|o Annealing : "annealing"
OneDMeasurement ||--|o Deposition : "deposition"
OneDMeasurement ||--|o Instrument : "instrument"
OneDMeasurement ||--|o Sample : "sample"
OneDMeasurement ||--|o TipSensor : "tip_sensor"
OneDMeasurement ||--|o TwoDMeasurement : "two_d_measurement"
TwoDMeasurement ||--|o QuantityValue : "duration"
TwoDMeasurement ||--|o QuantityValue : "bias_setpoint"
TwoDMeasurement ||--|o QuantityValue : "bias_calibration_factor"
TwoDMeasurement ||--|o QuantityValue : "bias_calibration_offset"
TwoDMeasurement ||--|o QuantityValue : "current_setpoint"
TwoDMeasurement ||--|o QuantityValue : "current_calibration_factor"
TwoDMeasurement ||--|o QuantityValue : "current_calibration_offset"
TwoDMeasurement ||--|o QuantityValue : "current_gain"
TwoDMeasurement ||--|o QuantityValue : "z_position"
TwoDMeasurement ||--|o QuantityValue : "z_controller_setpoint"
TwoDMeasurement ||--|o QuantityValue : "z_controller_time_constant"
TwoDMeasurement ||--|o QuantityValue : "z_controller_tip_lift"
TwoDMeasurement ||--|o QuantityValue : "z_controller_switch_off_delay"
TwoDMeasurement ||--|o PiezoConfigurationSettings : "piezo_configuration_settings"
TwoDMeasurement ||--|o ScanSettings : "scan_settings"
TwoDMeasurement ||--|o QuantityValue : "dwell_time"
TwoDMeasurement ||--|o QuantityValue : "sample_temperature"
TwoDMeasurement ||--|o QuantityValue : "recording_temperature"
TwoDMeasurement ||--|o OscillationControlSettings : "oscillation_control_settings"
TwoDMeasurement ||--|o QuantityValue : "e_min"
TwoDMeasurement ||--|o QuantityValue : "e_max"
TwoDMeasurement ||--|o QuantityValue : "de"
TwoDMeasurement ||--|o QuantityValue : "fwhm"
TwoDMeasurement ||--|o QuantityValue : "extrap_plane"
TwoDMeasurement ||--|o QuantityValue : "constant_height"
TwoDMeasurement ||--|o QuantityValue : "constant_current"
TwoDMeasurement ||--|o QuantityValue : "scan_dx"
TwoDMeasurement ||--|o QuantityValue : "scan_z_min"
TwoDMeasurement ||--|o QuantityValue : "scan_z_max"
TwoDMeasurement ||--|o QuantityValue : "amplitude"
TwoDMeasurement ||--|o QuantityValue : "resonance_frequency"
TwoDMeasurement ||--|o Annealing : "annealing"
TwoDMeasurement ||--|o Deposition : "deposition"
TwoDMeasurement ||--|o DFT : "dft"
TwoDMeasurement ||--|o TightBinding : "tb"
TwoDMeasurement ||--|o Instrument : "instrument"
TwoDMeasurement ||--|o Sample : "sample"
TwoDMeasurement ||--|o TipSensor : "tip_sensor"
Instrument ||--}o UHVComponent : "uhv_components"
Instrument ||--|o Room : "room"
Instrument ||--|o Manufacturer : "manufacturer"
Instrument ||--|o TipSensor : "tip_sensor"
DFT ||--|o QuantityValue : "scf_convergence_threshold"
DFT ||--|o QuantityValue : "absolute_magnetisation"
DFT ||--|o QuantityValue : "total_magnetisation"
DFT ||--|o QuantityValue : "smearing"
DFT ||--|o QuantityValue : "fermi_dirac_temperature"
DFT ||--|o OutputPopulationAnalysis : "output_mulliken_population_analysis"
DFT ||--|o OutputPopulationAnalysis : "output_hirshfeld_population_analysis"
Deposition ||--|o QuantityValue : "stabilisation_time"
Deposition ||--|o QuantityValue : "deposition_time"
Deposition ||--|o QuantityValue : "pressure"
Deposition ||--|o QuantityValue : "substrate_temperature"
Deposition ||--|o QuantityValue : "molecule_temperature"
Deposition ||--|o EvaporatorSlot : "evaporator_slot"
Deposition ||--|o Annealing : "annealing"
Deposition ||--|o Molecule : "molecule"
Deposition ||--}o UHVComponent : "uhv_components"
Deposition ||--|o Instrument : "instrument"
Deposition ||--|o Sample : "sample"
Molecule ||--}o EvaporationTemperature : "evaporation_temperatures"
Molecule ||--|o QuantityValue : "amount"
Molecule ||--|o Chemist : "chemist"
Molecule ||--|o Storage : "storage"
Molecule ||--}o Molecule : "precursor_molecules"
EvaporationTemperature ||--|o QuantityValue : "temperature"
Annealing ||--|o QuantityValue : "duration"
Annealing ||--|o QuantityValue : "pressure"
Annealing ||--|o QuantityValue : "voltage"
Annealing ||--|o QuantityValue : "temperature"
Annealing ||--|o QuantityValue : "current"
Annealing ||--|o Sputtering : "sputtering"
Annealing ||--|o Deposition : "deposition"
Annealing ||--}o UHVComponent : "uhv_components"
Annealing ||--|o Instrument : "instrument"
Annealing ||--|o Sample : "sample"
Sputtering ||--|o QuantityValue : "duration"
Sputtering ||--|o QuantityValue : "pressure"
Sputtering ||--|o QuantityValue : "discharge_voltage"
Sputtering ||--|o QuantityValue : "voltage"
Sputtering ||--|o QuantityValue : "temperature"
Sputtering ||--|o QuantityValue : "angle"
Sputtering ||--|o QuantityValue : "current"
Sputtering ||--}o UHVComponent : "uhv_components"
Sputtering ||--|o Annealing : "annealing"
Sputtering ||--|o Deposition : "deposition"
Sputtering ||--|o Instrument : "instrument"
Sputtering ||--|o Sample : "sample"
OscillationControlSettings ||--|o QuantityValue : "input_calibration"
OscillationControlSettings ||--|o QuantityValue : "input_range"
OscillationControlSettings ||--|o QuantityValue : "center_frequency"
OscillationControlSettings ||--|o QuantityValue : "range"
OscillationControlSettings ||--|o DemodulationSettings : "demodulation_settings"
OscillationControlSettings ||--|o PhaseSettings : "phase_settings"
OscillationControlSettings ||--|o QuantityValue : "frequency_shift"
OscillationControlSettings ||--|o AmplitudeSettings : "amplitude_settings"
OscillationControlSettings ||--|o QuantityValue : "excitation"
OscillationControlSettings ||--|o OutputSettings : "output_settings"
OscillationControlSettings ||--|o PLLSetupSettings : "pll_setup_settings"
PLLSetupSettings ||--|o QuantityValue : "demodulation_bw_amp"
PLLSetupSettings ||--|o QuantityValue : "demodulation_bw_pha"
PLLSetupSettings ||--|o QuantityValue : "amplitude_excitation_ratio"
OutputSettings ||--|o QuantityValue : "amplitude_range"
AmplitudeSettings ||--|o QuantityValue : "setpoint"
AmplitudeSettings ||--|o QuantityValue : "p_gain"
AmplitudeSettings ||--|o QuantityValue : "i_gain"
PhaseSettings ||--|o QuantityValue : "p_gain"
PhaseSettings ||--|o QuantityValue : "i_gain"
DemodulationSettings ||--}o QuantityValue : "frequencies"
DemodulationSettings ||--}o QuantityValue : "reference_phases"
DemodulationSettings ||--}o QuantityValue : "cutoff_frequencies"
DemodulationSettings ||--}o QuantityValue : "harmonics"
ScanSettings ||--|o QuantityValue : "angle"
ScanSettings ||--|o ScanOffset : "scan_offset"
ScanSettings ||--|o ScanTime : "scan_time"
ScanSettings ||--|o ScanRange : "scan_range"
ScanSettings ||--|o ScanPixels : "scan_pixels"
ScanSettings ||--|o ScanSpeed : "scan_speed"
ScanSpeed ||--|o QuantityValue : "forward_speed"
ScanSpeed ||--|o QuantityValue : "backward_speed"
ScanPixels ||--|o QuantityValue : "x_coordinate"
ScanPixels ||--|o QuantityValue : "y_coordinate"
ScanRange ||--|o QuantityValue : "x_coordinate"
ScanRange ||--|o QuantityValue : "y_coordinate"
ScanTime ||--|o QuantityValue : "x_coordinate"
ScanTime ||--|o QuantityValue : "y_coordinate"
ScanOffset ||--|o QuantityValue : "x_coordinate"
ScanOffset ||--|o QuantityValue : "y_coordinate"
PiezoConfigurationSettings ||--|o SensitivitySettings : "sensitivity_settings"
PiezoConfigurationSettings ||--|o HVGainSettings : "hv_gain_settings"
PiezoConfigurationSettings ||--|o ScanSlope : "scan_slope"
PiezoConfigurationSettings ||--|o CurvatureRadius : "curvature_radius"
PiezoConfigurationSettings ||--|o SecondOrderCorrection : "second_order_correction"
PiezoConfigurationSettings ||--|o DriftSettings : "drift_settings"
DriftSettings ||--|o QuantityValue : "x_coordinate"
DriftSettings ||--|o QuantityValue : "y_coordinate"
DriftSettings ||--|o QuantityValue : "z_coordinate"
SecondOrderCorrection ||--|o QuantityValue : "x_coordinate"
SecondOrderCorrection ||--|o QuantityValue : "y_coordinate"
CurvatureRadius ||--|o QuantityValue : "x_coordinate"
CurvatureRadius ||--|o QuantityValue : "y_coordinate"
ScanSlope ||--|o QuantityValue : "x_coordinate"
ScanSlope ||--|o QuantityValue : "y_coordinate"
HVGainSettings ||--|o QuantityValue : "x_coordinate"
HVGainSettings ||--|o QuantityValue : "y_coordinate"
HVGainSettings ||--|o QuantityValue : "z_coordinate"
SensitivitySettings ||--|o QuantityValue : "x_coordinate"
SensitivitySettings ||--|o QuantityValue : "y_coordinate"
SensitivitySettings ||--|o QuantityValue : "z_coordinate"
SRD ||--|o Molecule : "molecule"
Result ||--}o OneDMeasurement : "one_d_measurements"
Result ||--}o TwoDMeasurement : "two_d_measurements"
Result ||--}o GeometryOptimisation : "geometry_optimisations"
Result ||--}o Software : "softwares"
GeometryOptimisation ||--|o ForceConvergenceThreshold : "force_convergence_threshold"
GeometryOptimisation ||--}o AtomisticModel : "atomistic_models"
GeometryOptimisation ||--|o DFT : "dft"
GeometryOptimisation ||--|o EmpiricalModelling : "empirical_modelling"
AtomisticModel ||--}o AtomsPositions : "atoms_positions"
AtomisticModel ||--}o QuantityValue : "cell_vectors"
AtomisticModel ||--|o PeriodicBoundaryConditions : "periodic_boundary_conditions"
AtomisticModel ||--|o Crystal : "crystal"
AtomisticModel ||--|o Molecule : "molecule"
PeriodicBoundaryConditions ||--|o QuantityValue : "x_coordinate"
PeriodicBoundaryConditions ||--|o QuantityValue : "y_coordinate"
PeriodicBoundaryConditions ||--|o QuantityValue : "z_coordinate"
Publication ||--}o Grant : "grants"
Publication ||--}o Draft : "drafts"
Publication ||--}o Author : "authors"
Draft ||--}o Result : "results"
Grant ||--|o QuantityValue : "budget"
PDOS ||--}o QuantityValue : "energies"
PDOS ||--|o QuantityValue : "amplitude"
PDOS ||--|o QuantityValue : "degauss_energy"
PDOS ||--}o AtomisticModel : "atomistic_models"
PDOS ||--|o DFT : "dft"
PDOS ||--|o TightBinding : "tb"
MinimumEnergyPotential ||--|o QuantityValue : "energy_barrier"
MinimumEnergyPotential ||--}o QuantityValue : "energies"
MinimumEnergyPotential ||--}o AtomisticModel : "atomistic_models"
Maintenance ||--}o UHVComponent : "uhv_components"
Maintenance ||--|o Instrument : "instrument"
HydrogenCracker ||--}o UHVComponent : "uhv_components"
HydrogenCracker ||--|o Instrument : "instrument"
FillCryostat ||--|o QuantityValue : "weight_before"
FillCryostat ||--|o QuantityValue : "weight_after"
FillCryostat ||--}o UHVComponent : "uhv_components"
FillCryostat ||--|o Instrument : "instrument"
Errors ||--|o Instrument : "instrument"
Errors ||--}o UHVComponent : "uhv_components"
Dosing ||--|o QuantityValue : "duration"
Dosing ||--|o QuantityValue : "pressure"
Dosing ||--|o QuantityValue : "temperature"
Dosing ||--|o Instrument : "instrument"
Dosing ||--|o Sample : "sample"
CanGasBottle ||--|o QuantityValue : "pressure"
CanGasBottle ||--|o Supplier : "supplier"
BandStructure ||--|o KPointsConditions : "k_points_conditions"
BandStructure ||--|o BSEnergy : "bs_energies"
BandStructure ||--}o AtomisticModel : "atomistic_models"
BandStructure ||--|o DFT : "dft"
BandStructure ||--|o TightBinding : "tb"

```

