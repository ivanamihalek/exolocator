#!/usr/bin/env python

import svgwrite
import sys
from typing import Dict

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from typing_extensions import List, Dict

from el_utils.mysql import mysql_using_env_creds, mysql_server_conn_close
from el_utils.tree import species_tree, Node, Tree
from el_utils.bioinfo import get_seq_record_by_id

# IBM color blind palette
# https://www.color-hex.com/color-palette/1044488
'''
yello   #ffb000 	(255,176,0)
orange	#fe6100 	(254,97,0)
magenta	#dc267f 	(220,38,127)
purplish #785ef0 	(120,94,240)
powder  #648fff 	(100,143,255)

'''
tax_grouop2color = {
    'Primates':         '#648fff',
    'Euarchontoglires': '#785ef0',
    'Mammalia':         '#dc267f',
    'Tetrapoda':        '#fe6100',
    'Vertebrates':       '#ffb000'
}

legend_color_map = {
    'All vertebrates':      '#ffb000',
    'Tetrapods':            '#fe6100',
    'Mammals':              '#dc267f',
    'Rodents and primates': '#785ef0',
    'Primates':             '#648fff',
    'Variable in primates': '#ffffff',
}

def find_majority_character(characters: list[str]) -> str|None:
    """
    Finds the most frequent character(s) in an alignment column.

    Args:
        alignment_column: A string representing a single column from the alignment.

    Returns:
        A list of characters that appear most frequently in the column.
        If there's a tie, all characters with the maximum count are returned.
    """
    char_counts: Dict[str, int] = {}
    for char in characters:
        char_counts[char] = char_counts.get(char, 0) + 1
    max_count = max(char_counts.values())

    top_characters: List[str] = [c for c, count in char_counts.items() if count == max_count]
    if len(top_characters) == 1 and max_count/len(characters) >= 0.9:
        return top_characters[0]

    return None


def color_code(alignment: MultipleSeqAlignment, species_in_subtree: Dict[str, set]):
    # for each column i in the alignment
    color = {}
    human_sequence = str(get_seq_record_by_id(alignment, 'homo_sapiens').seq)
    ret_list = []
    for i in range(alignment.get_alignment_length()):
        color[i] = "#ffffff"
        # there is no guarantee how many seqs will be in each taxonomy group,
        # so I cannot just test is there is 90% conservation across all verterbrates
        # because it could be 100% conservation in Tetrapods, and none in the rest
        prev_set = set()
        for taxonomy_group, species in species_in_subtree.items():
            selected_chars = [record.seq[i] for record in alignment if record.id in species.difference(prev_set)]
            prev_set = species
            majority_character = find_majority_character(selected_chars)
            if not majority_character: break
            color[i] = tax_grouop2color[taxonomy_group]

        ret_list.append([i+1, human_sequence[i],  color[i]])

    return ret_list


def create_colored_stripe_svg(data: List[List], output_filename="colored_stripe.svg", font_size=14, blocks_per_line=50):
    """
    Generates an SVG file representing a colored stripe based on the input data.

    Args:
        input_data (str): A string containing lines of data in the format:
                           "index character color_hex".
                           Example: "1 M #ffb000"
        output_filename (str): The name of the SVG file to be created.
        font_size (int): The font size for the characters and numbers.
        blocks_per_line (int): Number of color blocks to display per line in the stripe.
    """

    padding_above = font_size
    block_height  = int(1.5 * font_size)  # Stripe height is 3 times the font size
    block_width   = 15  # Width of each color block
    x_offset      = 30  # horizontal offset
    line_spacing   = font_size  # Space between the stripe and the letters
    number_spacing = font_size  # Space between the letters and the index

    image_width = blocks_per_line * block_width + 2 * x_offset
    num_lines = (len(data) + blocks_per_line - 1) // blocks_per_line
    image_height = num_lines * (block_height + line_spacing + font_size + number_spacing + padding_above)  + padding_above

    svg_content = f'''<svg width="{image_width}" height="{image_height}" xmlns="http://www.w3.org/2000/svg">
    <defs>
        <style type="text/css">
            .text {{
                font-size: {font_size}px;
                font-family: monospace;
                font-weight: bold;
            }}
        </style>
    </defs>
    <rect width="100%" height="100%" fill="white"/>
    '''

    for line_num in range(num_lines):
        start_index = line_num * blocks_per_line
        end_index = min((line_num + 1) * blocks_per_line, len(data))

        # Calculate the Y offset for the current line
        y_offset = padding_above + line_num * (block_height + line_spacing + font_size + number_spacing + padding_above)

        # Add colored blocks for the current line
        for i in range(start_index, end_index):
            index, char, color = data[i]
            x = (i - start_index) * block_width + x_offset
            svg_content += f'''<rect x="{x}" y="{y_offset}" width="{block_width}" height="{block_height}" fill="{color}" />'''

            # Add character below the stripe
            text_x = x + block_width / 2
            text_y = y_offset + block_height + line_spacing + font_size / 2
            svg_content += f'''<text x="{text_x}" y="{text_y}" class="text" text-anchor="middle" fill="black">{char}</text>'''

            # Add index every tenth block
            if (index % 10 == 0):
                number_x = x + block_width / 2
                number_y = y_offset + block_height + line_spacing + font_size + number_spacing
                svg_content += f'''<text x="{number_x}" y="{number_y}" class="text" text-anchor="middle" fill="grey">{index}</text>'''

    svg_content += '</svg>'

    with open(output_filename, "w") as f:
        f.write(svg_content)

    print(f"SVG file created: {output_filename}")


def create_legend_svg(color_map, filename="legend.svg", width=400, height=200, box_size=20, margin=20,
                      border_color='gray', border_width=1):
    """
    Generates an SVG file containing a legend based on the provided color map.

    Args:
        color_map (dict): A dictionary mapping labels (strings) to colors (strings).
        filename (str): The name of the SVG file to create.
        width (int): The width of the SVG canvas.
        height (int): The height of the SVG canvas.
        box_size (int): The size (side length) of the color boxes in the legend.
        margin (int): Margin around the legend.
        border_color (str): The color of the rectangle borders.
        border_width (int): The width of the rectangle borders.

    """
    dwg = svgwrite.Drawing(filename, size=(width, height))
    x = margin
    y = margin
    line_height = box_size * 1.5  # Adjust line height based on box size

    for label, color in color_map.items():
        # Draw the color box
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=color,
                         stroke=border_color, stroke_width=border_width))

        # Add the label
        dwg.add(dwg.text(label, insert=(x + box_size + 10, y + box_size/1.3), fill='black', font_size=box_size/1.3))

        y += line_height  # Move to the next line

        if y + margin > height:
            print("Warning: Legend may exceed the specified height.")
            break

    dwg.save()
    print(f"Legend saved to {filename}")


def main():
    """
    Main function to read alignment, restricts the alignment to the target sequence,
    remove gappy sequences, and  write the output.
    """
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file> <output_file> ")
        sys.exit(1)

    [input_file, output_fnm] = sys.argv[1:3]

    try:
        # Read the alignment from the input file
        alignment = AlignIO.read(input_file, "fasta")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    species_list = [record.id for record in alignment]

    cursor = mysql_using_env_creds()
    tree: Tree = species_tree(cursor, species_list)
    mysql_server_conn_close(cursor)

    special_node_names = ['Primates', 'Euarchontoglires', 'Mammalia', 'Tetrapoda']
    species_in_subtree = {}
    for node_name in special_node_names:
        if not (parent_node := tree.get_node(node_name)): continue
        species_in_subtree[node_name] = set(parent_node.subtree_leaf_names())
    if len(species_list) > len(species_in_subtree['Tetrapoda']):
        species_in_subtree['Vertebrates'] = set(species_list)

    # assign color code to each position in the alignment according to conservation
    color_rep = color_code(alignment, species_in_subtree)
    create_colored_stripe_svg(color_rep, output_filename=output_fnm)
    create_legend_svg(legend_color_map)


###################################
if __name__ == "__main__":
    main()
