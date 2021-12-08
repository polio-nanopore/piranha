
# DEPENDENCIES AND RESOURCES TO CHECK
dependency_list = ["gofasta","minimap2","snakemake","iqtree","jclusterfunk"]
module_list = ["mako","Bio"]

resources = [
        {RESOURCE_KEY:"reference_sequence",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"references.fasta"},
        {RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data/report_modules",
        RESOURCE_KEY_FILENAME:"report_template.mako"}
    ]