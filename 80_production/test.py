from typing import List


def create_colored_stripe_svg(data: List[List], output_filename="colored_stripe.svg", font_size=12, blocks_per_line=50):
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

    padding_above = 2 * font_size
    block_height  = 3 * font_size  # Stripe height is 3 times the font size
    block_width   = 20  # Width of each color block
    x_offset = 30  # horizontal offset
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



# Example Usage (assuming your input data is in a string variable called 'input_data'):
input_data = """
1 M #ffb000
2 A #fe6100
3 S #785ef0
4 S #fe6100
5 P #fe6100
6 L #dc267f
7 P #fe6100
8 G #dc267f
9 P #dc267f
10 N #dc267f
11 D #dc267f
12 I #dc267f
13 L #fe6100
14 L #dc267f
15 A #dc267f
16 S #dc267f
17 P #dc267f
18 S #dc267f
19 S #dc267f
20 A #dc267f
21 F #785ef0
22 Q #dc267f
23 P #dc267f
24 D #785ef0
25 T #ffffff
26 L #dc267f
27 S #dc267f
28 Q #dc267f
29 P #648fff
30 R #dc267f
31 P #648fff
32 G #dc267f
33 H #dc267f
34 A #dc267f
35 N #785ef0
36 L #785ef0
37 K #dc267f
38 P #ffb000
39 N #ffb000
40 Q #ffb000
41 V #ffb000
42 G #dc267f
43 Q #ffb000
44 V #ffb000
45 I #fe6100
46 L #ffb000
47 Y #ffb000
48 G #ffb000
49 I #ffb000
50 P #ffb000
51 I #ffb000
52 V #ffb000
53 S #ffb000
54 L #ffb000
55 V #ffb000
56 I #ffb000
57 D #ffb000
58 G #fe6100
59 Q #fe6100
60 E #ffb000
61 R #ffb000
62 L #ffb000
63 C #ffb000
64 L #ffb000
65 A #ffb000
66 Q #ffb000
67 I #ffb000
68 S #ffb000
69 N #ffb000
70 T #ffb000
71 L #ffb000
72 L #ffb000
73 K #ffb000
74 N #ffb000
75 F #dc267f
76 S #ffb000
77 Y #ffb000
78 N #ffb000
79 E #ffb000
80 I #ffb000
81 H #ffb000
82 N #ffb000
83 R #ffb000
84 R #ffb000
85 V #ffb000
86 A #ffb000
87 L #ffb000
88 G #ffb000
89 I #ffb000
90 T #ffb000
91 C #ffb000
92 V #ffb000
93 Q #ffb000
94 C #ffb000
95 T #ffb000
96 P #ffb000
97 V #ffb000
98 Q #ffb000
99 L #ffb000
100 E #ffb000
101 I #ffb000
102 L #ffb000
103 R #ffb000
104 R #ffb000
105 A #ffb000
106 G #ffb000
107 A #ffb000
108 M #ffb000
109 P #ffb000
110 I #ffb000
111 S #ffb000
112 S #ffb000
113 R #ffb000
114 R #ffb000
115 C #ffb000
"""

lines = input_data.strip().split('\n')
data = []
for line in lines:
    parts = line.split()
    if len(parts) == 3:
        index, char, color = parts
        data.append((int(index), char, color))

create_colored_stripe_svg(data)
