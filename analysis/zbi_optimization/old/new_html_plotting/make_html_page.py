import os
import argparse
import glob as glob

def int_sorting_key(filename):
    parts = filename.split('_')
    if parts[1].isdigit():
        return (int(parts[1]), filename)
    else:
        return (1, filename)  # For non-integer keys, keep them in their original order

# Function to extract a key for sorting
def sorting_key(filename):
    parts = filename.split('_')
    if parts[1] == "initial":
        return (0, filename)
    elif parts[1].isdigit():
        return (int(parts[1]), filename)
    else:
        return (1, filename)  # For non-integer keys, keep them in their original order

def make_parent_html(base_dir, name):
    print("making dir", base_dir)
    os.makedirs(base_dir, exist_ok=True)
    with open(os.path.join(base_dir, '%s_parent.html'%(name)), 'w') as f:
        f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{name} ZBi Processor Results</title>
    </head>
    <body>
        <h1> {name} ZBi Plots</h1>
        <p>Click on the links below to explore variables and iterations:</p>
        
        <ul>
            <li><a href="variables.html">variables</a></li>
            <li><a href="iterations.html">iterations</a></li>
        </ul>
    </body>
    </html>
        ''')

def make_master_html(base_dir, name):
    pages = [page for page in sorted(os.listdir(base_dir))]
    html_files = glob.glob(os.path.join(base_dir, "*/*parent.html"))
    html_files = [file.replace(base_dir+'/','') for file in html_files]
    print("LOOOOOOK: ", html_files)
    print("base_dir:" , base_dir)

    #Write html file with a link to the html page of each variable
    with open(os.path.join(base_dir, '%s.html'%(name)), 'w') as f:
        #Automatically generate links for every variable in list of variables
        links = "\n".join([f'<li><a href="{html}">{html}</a></li>' for html in html_files])
        f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{top_name} - Variables - 1D Plots As Function of Iteration</title>
    </head>
    <body>
        <h1>Signal and Background Variables for {top_name}</h1>
        <p> Click a variable to see how Signal and Background change with iterative cuts </p>
        <ul>
            {links}
        </ul>
    </body>
    </html>
        ''')

def make_iterations_html(parent_dir, top_name, name):
    base_dir = os.path.join(parent_dir, name)
    basename = os.path.basename(base_dir)
    os.makedirs(base_dir, exist_ok=True)
    with open(os.path.join(parent_dir, '%s.html'%(name)), 'w') as f:
        f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{top_name} Page</title>
    </head>
    <body>
        <h1>{top_name} - {name} </h1>
        <p>Click:</p>
        
        <ul>
            <li><a href="{basename}/background_models.html">1D Plots</a></li>
        </ul>
    </body>
    </html>
        ''')

    return base_dir

def make_variables_html(parent_dir, top_name, name):
    base_dir = os.path.join(parent_dir, name)
    basename = os.path.basename(base_dir)
    os.makedirs(base_dir, exist_ok=True)
    with open(os.path.join(parent_dir, '%s.html'%(name)), 'w') as f:
        f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{top_name} Page</title>
    </head>
    <body>
        <h1>{top_name} - {name} </h1>
        <p>Click 1D or 2D plots:</p>
        
        <ul>
            <li><a href="{basename}/plots_1D.html">1D Plots</a></li>
            <li><a href="{basename}/plots_2D.html">2D Plots</a></li>
        </ul>
    </body>
    </html>
        ''')

    return base_dir

def make_var_1d_plots(parent_dir, top_name, name, input_images_dir, keyword):

    #Get list of available variables from image names
    #input_images_dir = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/zbi_html/zalpha_plots_html_test/1d_plots'
    variables = [png for png in sorted(os.listdir(input_images_dir)) if keyword in png]
    variables = [os.path.basename(png).replace(keyword,'').replace('.png','') for png in variables]
    print("List of available variables: ", variables)

    #Make 1D Plot subdirectory
    base_dir = os.path.join(parent_dir, name)
    basename = os.path.basename(base_dir)
    os.makedirs(base_dir, exist_ok=True)

    #Write html file with a link to the html page of each variable
    with open(os.path.join(parent_dir, '%s.html'%(name)), 'w') as f:
        #Automatically generate links for every variable in list of variables
        links = "\n".join([f'<li><a href="{basename}/{var}.html">{var}</a></li>' for var in variables])
        f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{top_name} - Variables - 1D Plots As Function of Iteration</title>
    </head>
    <body>
        <h1>Signal and Background Variables for {top_name}</h1>
        <p> Click a variable to see how Signal and Background change with iterative cuts </p>
        <ul>
            {links}
        </ul>
    </body>
    </html>
        ''')

    #Write an html file for each variable that contains links to variable plots
    for var in variables:
        with open(os.path.join(base_dir, f'{var}.html'), 'w') as f:

            #Generate links to variable plots
            matching_pngs = [png for png in sorted(os.listdir(input_images_dir)) if var in png]
            sorted_pngs = sorted(matching_pngs, key=sorting_key)
            # Generate div containers with clickable img tags for the sorted PNG files
            img_divs = "\n".join([
                f'<div><a href="{input_images_dir}/{png}" target="_blank"><img src="{input_images_dir}/{png}" alt="{png}" width="1000"></a></div>'
                for png in sorted_pngs
            ])

            f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{top_name} {var} 1D Plots By Iteration</title>
        <style>
            div {{
                margin-bottom: 20px; /* Add spacing between images */
            }}
        </style>
    </head>
    <body>
        <h1></h1>
        <ul>
            {img_divs}
        </ul>
    </body>
    </html>
            ''')


def make_var_2d_plots(parent_dir, top_name, name, input_images_dir, keyword):

    #Get list of available variables from image names
    variables = [png for png in sorted(os.listdir(input_images_dir)) if keyword in png]
    variables = [os.path.basename(png).replace(keyword,'').replace('.png','') for png in variables]
    variables.append("recon_z")
    print("List of available variables: ", variables)

    #Make Plot subdirectory
    base_dir = os.path.join(parent_dir, name)
    basename = os.path.basename(base_dir)
    os.makedirs(base_dir, exist_ok=True)

    #Write html file with a link to the html page of each variable
    with open(os.path.join(parent_dir, '%s.html'%(name)), 'w') as f:
        #Automatically generate links for every variable in list of variables
        links = "\n".join([f'<li><a href="{basename}/{var}.html">{var}</a></li>' for var in variables])
        f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{top_name} - Variables - 2D Plots As Function of Iteration</title>
    </head>
    <body>
        <h1>Signal and Background Variables for {top_name}</h1>
        <p> Click a variable to see how Signal and Background change with iterative cuts </p>
        <ul>
            {links}
        </ul>
    </body>
    </html>
        ''')

    #Write an html file for each variable that contains links to variable plots
    for var in variables:
        with open(os.path.join(base_dir, f'{var}.html'), 'w') as f:

            #Generate links to variable plots
            matching_pngs = [png for png in sorted(os.listdir(input_images_dir)) if var in png]
            sorted_pngs = sorted(matching_pngs, key=sorting_key)
            # Generate div containers with clickable img tags for the sorted PNG files
            img_divs = "\n".join([
                f'<div><a href="{input_images_dir}/{png}" target="_blank"><img src="{input_images_dir}/{png}" alt="{png}" width="1000"></a></div>'
                for png in sorted_pngs
            ])

            f.write(f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{top_name} {var} 2D Plots By Iteration</title>
        <style>
            div {{
                margin-bottom: 20px; /* Add spacing between images */
            }}
        </style>
    </head>
    <body>
        <h1></h1>
        <ul>
            {img_divs}
        </ul>
    </body>
    </html>
            ''')

def make_bkg_model_plots(parent_dir, top_name, name, input_images_dir):

    #Generate links to variable plots
    sorted_pngs = [png for png in sorted(os.listdir(input_images_dir))]
    sorted_pngs = sorted(sorted_pngs, key=int_sorting_key)

    #Make Plot subdirectory
    base_dir = os.path.join(parent_dir, name)
    basename = os.path.basename(base_dir)
    os.makedirs(base_dir, exist_ok=True)

    #Write html file with a link to the html page of each variable
    with open(os.path.join(parent_dir, '%s.html'%(name)), 'w') as f:
        # Generate div containers with clickable img tags for the sorted PNG files
        img_divs = "\n".join([
            f'<div><a href="{input_images_dir}/{png}" target="_blank"><img src="{input_images_dir}/{png}" alt="{png}" width="1000"></a></div>'
            for png in sorted_pngs
        ])

        f.write(f'''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{top_name} Iteration Background Models</title>
    <style>
        div {{
            margin-bottom: 20px; /* Add spacing between images */
        }}
    </style>
</head>
<body>
    <h1></h1>
    <ul>
        {img_divs}
    </ul>
</body>
</html>
        ''')

###############################################################################################
parser = argparse.ArgumentParser(description="HTML Config")
parser.add_argument('--name', type=str, dest="name", default = "variable_html")
parser.add_argument('--dir', type=str, dest="base_dir", default = "html_test")
options = parser.parse_args()

# Directory path
plots_dir = "./plots_1d"

master_dir = "master_html_pages"

# List all directories that will become html pages
html_page_names = [d for d in os.listdir(plots_dir) if os.path.isdir(os.path.join(plots_dir, d))]

input_plots_1d_dir = '/home/alic/scpTrash/html_pages/plots_1d/'
input_plots_2d_dir = '/home/alic/scpTrash/html_pages/plots_2d/'
input_bkg_models_dir = "/home/alic/scpTrash/html_pages/bkg_models/"

# Loop over each directory
for html_page_name in html_page_names:
    print("Making html page for ", html_page_name)
    html_page_name = os.path.basename(html_page_name)
    top_name = html_page_name
    top_dir = master_dir+'/'+html_page_name

    #Parent
    make_parent_html(top_dir, top_name)

    #Variables
    var_dir = make_variables_html(top_dir, top_name, 'variables')

    #Variable 1D Plots
    #input_images_dir = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/zbi_html/zalpha_plots_html_test/1d_plots'
    make_var_1d_plots(var_dir, top_name, 'plots_1D', input_plots_1d_dir+top_name, "iteration_initial_")

    #Variable 2D Plots
    #input_images_dir = "/sdf/group/hps/users/alspellm/projects/THESIS/ana/zbi_html/zalpha_plots_html_test/2d_plots"
    make_var_2d_plots(var_dir, top_name, 'plots_2D', input_plots_2d_dir+top_name, "initial_")

    #Iteations
    iter_dir = make_iterations_html(top_dir, top_name, 'iterations')

    #Background Model Plots
    #input_images_dir = "/sdf/group/hps/users/alspellm/projects/THESIS/ana/zbi_html/zalpha_plots_html_test/2d_plots"
    make_bkg_model_plots(iter_dir, top_name, 'background_models', input_bkg_models_dir+top_name)

make_master_html(master_dir, 'master')
print("File structure and HTML pages have been created.")
