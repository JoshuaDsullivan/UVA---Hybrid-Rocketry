import numpy as np

def calculate_o_ring_parameters(od, id, thickness):
    groove_depths = np.linspace(.05, .1, .15, .5, 1, 100)  # Range of possible groove depths (adjust as needed)
    valid_groove_depths = []
    
    for groove_depth in groove_depths:
        # Step 1: Compute Percent Stretch
        initial_id = id
        stretched_id = groove_depth - 4.75  # Oxidizer tank inner diameter is fixed at 4.75"
        stretch_percent = (stretched_id - initial_id) / stretched_id * 100
        
        if not (1 <= stretch_percent <= 5):
            continue  # Skip groove depths that don't satisfy stretch condition
        
        # Step 2: Compute Percent Reduction in Cross Section
        Y = 100 ** (1 - 10 / np.sqrt(100 + stretch_percent))
        new_thickness = thickness * (1 - Y / 100)
        
        # Step 3: Compute New Inner Diameter
        new_inner_diameter = stretched_id * (1 - Y / 100)
        
        # Step 4: Compute New Outer Diameter
        new_outer_diameter = new_inner_diameter + 2 * new_thickness
        
        # Step 5: Compute Compression Ratio
        compression_ratio = (thickness - new_thickness) / thickness * 100
        
        if 8 <= compression_ratio <= 35:
            valid_groove_depths.append((groove_depth, stretch_percent, compression_ratio))
    
    return valid_groove_depths

# Example usage
o_ring_od = 4.728  # Outer diameter in inches
o_ring_id = 4.45  # Inner diameter in inches
o_ring_thickness = 0.139  # Thickness in inches

results = calculate_o_ring_parameters(o_ring_od, o_ring_id, o_ring_thickness)
for depth, stretch, compression in results:
    print(f"Groove Depth: {depth:.3f} in, Stretch: {stretch:.2f}%, Compression: {compression:.2f}%")
