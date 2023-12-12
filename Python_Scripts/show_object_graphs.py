import networkx as nx
import plotly.graph_objects as go

project_description = "First experiments for C60 on hBN on Rh(111)"

# Create a graph
G = nx.DiGraph()

# Add nodes
G.add_node("ANL1", 
           type="Annealing", name="524a on Au(111) anneal", description="", duration_minutes="15 min",
           temperature_celsius="410", current_amperes="1.1", voltage_volts="19.67", pressure_mbars="",
           comments="")

# G.add_node("AUTH1", 
#            type="Author", name="F. Lopes")

# G.add_node("BND_SIM1", 
#            type="Bands Simultation", name="", simulation_code="", comments="")

# G.add_node("UHV_COMP1", 
#            type="UHV Component", name="", uhv_component_type="", uhv_component_model="",
#            serial_number="", empa_id="", website="", receive_date="", comments="")

G.add_node("CRYS1", 
           type="Crystal", name="Au111-Okan", material="", face="",
           sample_plate="", diameter_mm="", height_mm="", 
           hazardous="", hazardous_specification="", fridge="",
           no_light="", dry="", no_ox="",
           other_storage_condition="", other_storage_condition_specification="", specifications="",
           receive_date="", reference_number="", comments="")

# G.add_node("DRF1", 
#            type="Draft", name="", draft_type="", 
#            comments="")

# G.add_node("EVAP1", 
#            type="Evaporation", name="", evaporation_type="", 
#            comments="")

# G.add_node("GRT1", 
#            type="Grant", name="", funding_agency="", 
#            start_date="", end_date="", budget="",
#            currency="", project_id="", acknowledgement_sentence="",
#            comments="")

G.add_node("INST1", 
           type="Institution", name="Empa", address="", 
           comments="")

G.add_node("INSTR1", 
           type="Instrument", name="SRD", serial_number="", 
           empa_id="", comments="")

G.add_node("INSTR2", 
           type="Instrument", name="QPlus FEL", serial_number="", 
           empa_id="", comments="")

# G.add_node("L2DM1", 
#            type="Layered 2D Material", name="", empa_number="",
#            batch="", sample_plate="", width_mm="", 
#            length_mm="", thickness_mm="", hazardous="", 
#            hazardous_specification="", fridge="", no_light="", 
#            dry="", no_ox="", other_storage_condition="", 
#            other_storage_condition_specification="", receive_date="", comments="")

G.add_node("LOC1", 
           type="Location", name="Air", location_type="", 
           comments="")

# G.add_node("MAINT1", 
#            type="Maintenance", name="", description="", 
#            comments="")

# G.add_node("MANUFPER1", 
#            type="Manufacturer Person", name="", email="", 
#            website="")

G.add_node("MOL1", 
           type="Molecule", name="524a", empa_number="",
           batch="", iupac_name="", sum_formula="", 
           smiles="", cas_number="", hazardous="", 
           hazardous_specification="", evaporation_temperature="", fridge="", 
           no_light="", dry="", no_ox="", 
           other_storage_condition="", other_storage_condition_specification="", molecule_vials="",
           molecule_crucibles="", amount="", receive_date="", 
           comments="")

# G.add_node("PROD1", 
#            type="Molecule Product", name="", empa_number="",
#            batch="", iupac_name="", sum_formula="", 
#            smiles="", cas_number="", comments="")

G.add_node("TRANSF1", 
           type="Transfer", name="Transfer of sample 524a", description="From SRD to Air", 
           comments="")

G.add_node("TRANSF2", 
           type="Transfer", name="Transfer of sample 524a", description="From Air to QPlus FEL", 
           comments="")

# G.add_node("OPTIC1", 
#            type="Optical Component", name="", optical_type="", 
#            optical_model_name="", serial_number="", website="",
#            receive_date="", comments="")

# G.add_node("PDOS1", 
#            type="PDOS Simulation", name="", simulation_code="", 
#            comments="")

G.add_node("PERS1", 
           type="Person", name="Amogh Kinikar", short_name="kiam", 
           work_status="Active", email="", work_phone="",
           mobile_phone="")

G.add_node("PREP1", 
           type="Preparation", name="524a on Au(111)", preparation_date="", 
           comments="")

# G.add_node("PUBL1", 
#            type="Publication", name="", abstract="", 
#            doi="", year="", oa_version="",
#            oa_data="", comments="")

# G.add_node("RES1", 
#            type="Result", name="", description="", 
#            comments="")

G.add_node("ROOM1", 
           type="Room", name="La148", comments="")

G.add_node("ROOM2", 
           type="Room", name="La042", comments="")

# G.add_node("SPMI1", 
#            type="SPM Image", name="", description="", 
#            comments="")

G.add_node("SPUT1", 
           type="Sputtering", name="524 on Au(111) sput", duration_minutes="10 min", 
           pressure_mbars="8E-6 mBar", energy_kilovolts="1 kV", voltage_volts="780 V", temperature_celsius="",
           angle_degrees="", current_milliamperes="6.8 mA", comments="")

# G.add_node("SUPPL1", 
#            type="Supplier", name="")

# G.add_node("TIP1", 
#            type="Tip Sensor", name="", tip_sensor_type="",
#            wire_material="", etched_fib="", holder_id="", 
#            thickness_mm="", length_mm="", resonance_frequency="", 
#            quality_factor="", reference_number="", specifications="", 
#            receive_date="", comments="")

# G.add_node("SUB1", 
#            type="Wafer Substrate", name="", empa_number="",
#            batch="", sample_number="", width_mm="", 
#            length_mm="", thickness_mm="", hazardous="", 
#            hazardous_specification="", fridge="", no_light="", 
#            dry="", no_ox="", other_storage_condition="", 
#            other_storage_condition_specification="", receive_date="", comments="")

# G.add_node("THDEV1", 
#            type="Theory Device", name="", description="", 
#            comments="")

# G.add_node("ERROR1", 
#            type="Error", name="", description="", 
#            comments="")

# G.add_node("ANLYS1", 
#            type="Analysis", name="", description="", 
#            comments="")

# G.add_node("CALCTION1", 
#            type="Calculation", name="", description="", 
#            comments="")

# G.add_node("THEORY_OTHER1", 
#            type="Theory Task Other", name="", description="", 
#            comments="")

# G.add_node("SLAB_CALC1", 
#            type="SLAB Calculation", name="", description="", 
#            comments="")

# G.add_node("REACT1", 
#            type="Reaction", name="", description="", 
#            simulation_code="", reaction_barrier="", comments="")

# G.add_node("UPLOAD_STRUCT1", 
#            type="Upload Structure", name="", description="", 
#            comments="")

# G.add_node("NANORIB_WORK1", 
#            type="Nanoribbon Workflow", name="", description="", 
#            comments="")

# G.add_node("CELL_OPT1", 
#            type="Cell Optimisation", name="", description="", 
#            comments="")

# G.add_node("MOL_GAS_PHS1", 
#            type="Molecule Gas Phase", name="", description="", 
#            comments="")

# G.add_node("CALIBTION1", 
#            type="Calibration Optimisation", name="", description="", 
#            comments="")

# G.add_node("FILL_CRYO1", 
#            type="Fill Cryostat", name="", description="", 
#            cryostat_liquid="", cryostat_weight_before="", cryostat_weight_after="",
#            comments="")

# G.add_node("SYST_STAT1", 
#            type="System Status", name="", description="", 
#            comments="")

# G.add_node("SYST_HNDOVR1", 
#            type="System Handover", name="", description="", 
#            comments="")

# G.add_node("DEGAS_MOL1", 
#            type="Degas Molecule", name="", description="", 
#            comments="")

# G.add_node("AFM1", 
#            type="AFM", name="", description="",
#            mode="", tip="", resonance_frequency="", 
#            scan_rate="", drive_amplitude="", proportional_gain="", 
#            integral_gain="", comments="")

# G.add_node("ARPES1", 
#            type="ARPES", name="", description="", 
#            comments="")

# G.add_node("KPFM1", 
#            type="KPFM", name="", description="",
#            scan_number="", voltage_range_min_volts="", voltage_range_max_volts="", 
#            num_points="", raster_time="", delay_times="", 
#            frequency_hertz="", oscillation_q="", amplitude_mv="", 
#            oscillation_amplitude_p="", oscillation_amplitude_i="", mx_pll_p="", 
#            mx_pll_i="", mx_pll_bw="", advisor_bw="",
#            advisor_max_peak="", phase_detector_bw="", pi_p="", 
#            pi_tau="", comments="")

# G.add_node("LEED1", 
#            type="LEED", name="", description="", 
#            comments="")

# G.add_node("MASS_SPEC1", 
#            type="Mass-Spec", name="", description="", 
#            comments="")

# G.add_node("NC_AFM1", 
#            type="Non-contact AFM", name="", description="",
#            scan_number="", frequency_hertz="", oscillation_q="", 
#            amplitude_mv="", oscillation_amplitude_p="", oscillation_amplitude_i="", 
#            mx_pll_p="", mx_pll_i="", mx_pll_bw="", advisor_bw="",
#            advisor_max_peak="", phase_detector_bw="", pi_p="", 
#            pi_tau="", comments="")

G.add_node("STM1", 
           type="STM", name="", description="Decent coverage. One species on surface",
           scan_number="", measurement_temperature="", comments="")

# G.add_node("STS1", 
#            type="STS", name="", description="",
#            scan_number="", frequency_hertz="", phase_degrees="", 
#            excitation_amplitude_mv="", calibration_factor_volts="", comments="")

# G.add_node("TDP1", 
#            type="TDP", name="", description="", 
#            comments="")

# G.add_node("UPS1", 
#            type="UPS", name="", description="",
#            source="", pressure_arp="", pressure_3rd="", 
#            pressure_2nd="", pressure_1st="", current_amperes="", 
#            voltage_volts="", angle_degrees="", peak_type="", peak_voltage_min="",
#            peak_voltage_max="", peak_voltage_step_size="", dwell_time="", 
#            pass_energy="", mode="", slit="", comments="")

# G.add_node("XPS1", 
#            type="XPS", name="", description="",
#            source="", pressure_arp="", pressure_3rd="", 
#            pressure_2nd="", pressure_1st="", current_amperes="", 
#            voltage_volts="", angle_degrees="", peak_type="", peak_voltage_min="",
#            peak_voltage_max="", peak_voltage_step_size="", dwell_time="", 
#            pass_energy="", mode="", slit="", comments="")

# G.add_node("XPD1", 
#            type="XPD", name="", description="",
#            source="", pressure_arp="", pressure_3rd="", 
#            pressure_2nd="", pressure_1st="", current_amperes="", 
#            voltage_volts="", angle_degrees="", peak_type="", peak_voltage_min="",
#            peak_voltage_max="", peak_voltage_step_size="", dwell_time="", 
#            pass_energy="", mode="", slit="", comments="")

G.add_node("QCM1", 
           type="QCM", name="QCM", description="",
           pressure_mbars="7e-10 mBar", rate_angstroem_min="0.3 Angstroem/min", temperature_celsius="215 C", 
           pid="40, 17, 3", comments="")

G.add_node("QCM2", 
           type="QCM", name="QCM", description="",
           pressure_mbars="7e-10 mBar", rate_angstroem_min="0.7 Angstroem/min", temperature_celsius="235 C", 
           pid="40, 17, 3", comments="")

# G.add_node("EXPER_OTHER1", 
#            type="Experiment Task Other", name="", description="", 
#            comments="")

# G.add_node("BAKEOUT1", 
#            type="Bakeout", name="", description="", 
#            comments="")

# G.add_node("SRD1", 
#            type="SRD", name="", description="", 
#            comments="")

G.add_node("DEPOS1", 
           type="Deposition", name="", description="", 
           duration_minutes="7 min+", deposition_pressure="7e-10 mBar", substrate_temperature="RT", 
           molecule_temperature="235 C", comments="")

# G.add_node("DOSING1", 
#            type="Dosing", name="", description="", 
#            duration_minutes="", pressure_mbars="", comments="")

# G.add_node("LIGHT_IRRAD1", 
#            type="Light irradiation", name="", description="", 
#            comments="")

# G.add_node("FIELD_EMIS1", 
#            type="Field emission", name="", description="", 
#            comments="")

# G.add_node("HYDRO_CRACK1", 
#            type="Hydrogen Cracker", name="", description="", 
#            comments="")

# G.add_node("CATALYS_CHAMB1", 
#            type="Catalysis Chamber", name="", description="", 
#            comments="")

# G.add_node("FIB1", 
#            type="FIB", name="", description="", 
#            comments="")

# G.add_node("SEM1", 
#            type="SEM", name="", description="", 
#            comments="")

# G.add_node("EDX1", 
#            type="EDX", name="", description="", 
#            comments="")

# G.add_node("RAMAN1", 
#            type="Raman", name="", description="", 
#            laser_wavelength="", raman_objective="", laser_power="", 
#            sweeps_integration_time="", area_length_um="", area_width_um="", 
#            polarisation="", vacuum="", comments="")

# G.add_node("UV_VIS1", 
#            type="UV-Vis", name="", description="", 
#            mode="", measurement_range_min="", measurement_range_max="", 
#            scan_speed="", comments="")

# G.add_node("NEB1", 
#            type="NEB", name="", description="", 
#            comments="")

# G.add_node("CHAIN_GEO_OPT1", 
#            type="Chain_geo_opt", name="", description="", 
#            comments="")

# Add edges to connect objects
G.add_edge("PERS1", "INST1")
G.add_edge("INSTR1", "ROOM1")
G.add_edge("INSTR2", "ROOM2")

G.add_edge("TRANSF1", "PERS1")
G.add_edge("TRANSF1", "LOC1")
G.add_edge("TRANSF1", "INSTR1")
G.add_edge("TRANSF1", "MOL1")

G.add_edge("TRANSF2", "PERS1")
G.add_edge("TRANSF2", "LOC1")
G.add_edge("TRANSF2", "INSTR2")
G.add_edge("TRANSF2", "MOL1")

G.add_edge("QCM1", "INSTR1")
G.add_edge("QCM2", "INSTR1")
G.add_edge("QCM1", "MOL1")
G.add_edge("QCM1", "QCM2")

G.add_edge("PREP1", "CRYS1")
G.add_edge("PREP1", "MOL1")
G.add_edge("PREP1", "INSTR1")
G.add_edge("PREP1", "SPUT1")
G.add_edge("SPUT1", "ANL1")
G.add_edge("ANL1", "DEPOS1")
G.add_edge("DEPOS1", "STM1")


# Define colors for different node types
color_names = [
    "aliceblue", "antiquewhite", "aqua", "aquamarine", "azure",
    "beige", "bisque", "black", "blanchedalmond", "blue",
    "blueviolet", "brown", "burlywood", "cadetblue",
    "chartreuse", "chocolate", "coral", "cornflowerblue",
    "cornsilk", "crimson", "cyan", "darkblue", "darkcyan",
    "darkgoldenrod", "darkgray", "darkgrey", "darkgreen",
    "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange",
    "darkorchid", "darkred", "darksalmon", "darkseagreen",
    "darkslateblue", "darkslategray", "darkslategrey",
    "darkturquoise", "darkviolet", "deeppink", "deepskyblue",
    "dimgray", "dimgrey", "dodgerblue", "firebrick",
    "floralwhite", "forestgreen", "fuchsia", "gainsboro",
    "ghostwhite", "gold", "goldenrod", "gray", "grey", "green",
    "greenyellow", "honeydew", "hotpink", "indianred", "indigo",
    "ivory", "khaki", "lavender", "lavenderblush", "lawngreen",
    "lemonchiffon", "lightblue", "lightcoral", "lightcyan",
    "lightgoldenrodyellow", "lightgray", "lightgrey",
    "lightgreen", "lightpink", "lightsalmon", "lightseagreen",
    "lightskyblue", "lightslategray", "lightslategrey",
    "lightsteelblue", "lightyellow", "lime", "limegreen",
    "linen", "magenta", "maroon", "mediumaquamarine",
    "mediumblue", "mediumorchid", "mediumpurple",
    "mediumseagreen", "mediumslateblue", "mediumspringgreen",
    "mediumturquoise", "mediumvioletred", "midnightblue",
    "mintcream", "mistyrose", "moccasin", "navajowhite", "navy",
    "oldlace", "olive", "olivedrab", "orange", "orangered",
    "orchid", "palegoldenrod", "palegreen", "paleturquoise",
    "palevioletred", "papayawhip", "peachpuff", "peru", "pink",
    "plum", "powderblue", "purple", "red", "rosybrown",
    "royalblue", "rebeccapurple", "saddlebrown", "salmon",
    "sandybrown", "seagreen", "seashell", "sienna", "silver",
    "skyblue", "slateblue", "slategray", "slategrey", "snow",
    "springgreen", "steelblue", "tan", "teal", "thistle", "tomato",
    "turquoise", "violet", "wheat", "white", "whitesmoke",
    "yellow", "yellowgreen"
]

node_colors = {G.nodes()[node]['type']: color_names[i] for i, node in enumerate(G.nodes())}

# Create a Plotly graph
pos = nx.spring_layout(G)  # You can use other layout options
edge_trace = go.Scatter(
    x=[],
    y=[],
    line=dict(width=0.5, color="#888"),
    hoverinfo="none",
    mode="lines")

node_trace = go.Scatter(
    x=[],
    y=[],
    text=[],
    mode="markers+text",
    hoverinfo="text",
    marker=dict(
        showscale=False,
        size=20,
        color=[node_colors[data["type"]] for node, data in G.nodes(data=True)],
        line_width=2))

hovertext = []
for node in G.nodes():
    x, y = pos[node]
    node_trace["x"] += (x,)
    node_trace["y"] += (y,)
    node_trace["text"] += (node,)

    node_attributes = G.nodes[node]
    hovertext.append(node_attributes)
    
node_trace["hovertext"] = hovertext

for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_trace["x"] += (x0, x1, None)
    edge_trace["y"] += (y0, y1, None)

fig = go.Figure(
    data=[edge_trace, node_trace],
    layout=go.Layout(
        showlegend=False,
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))

)

# Set node colors based on type
node_color_list = [node_colors[data["type"]] for node, data in G.nodes(data=True)]
node_trace.marker.color = node_color_list

# Add attributes to hover text
hover_text = []
for node, data in G.nodes(data=True):
    hover_text.append(f"{node}<br>Name: {data['name']}")

node_trace.text = hover_text

fig.show()
