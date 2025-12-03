# Import libraries
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import seaborn as sns
import os
from PIL import Image


def plot_celiac_samples_over_time(all_samples_df, output_dir, write_tsv_path=None):
    """
    Number of diagnosed celiac samples available over time.

    This version fills in every month and uses a solid step‐style area to
    illustrate dataset releases.
    """
    # Mapping for month names to standard abbreviations
    month_map = {
        'Sept': 'Sep',
        'Sept-': 'Sep-'
    }

    # Standardize month abbreviations
    all_samples_df['Month_of_Publication'] = all_samples_df['Month_of_Publication'].replace(month_map, regex=True)

    # Convert Month_of_Publication to datetime (first day of month)
    all_samples_df['Month_of_Publication'] = pd.to_datetime(all_samples_df['Month_of_Publication'], format='%b-%y')

    # Sort by date
    all_samples_df = all_samples_df.sort_values('Month_of_Publication')

    # Count cumulative celiac samples over time
    celiac_samples = (
        all_samples_df[all_samples_df['Diagnosed_Celiac'] == True]
        .groupby('Month_of_Publication')
        .size()
        .cumsum()
    )

    # Add zero point one month before the first date
    first_date = celiac_samples.index.min()
    zero_date = first_date - pd.DateOffset(months=1)
    celiac_samples = pd.concat([pd.Series({zero_date: 0}), celiac_samples])

    # Ensure every month is represented
    last_date = celiac_samples.index.max()
    full_months = pd.date_range(start=zero_date, end=last_date, freq='MS')
    celiac_samples = (
        celiac_samples
        .reindex(full_months)
        .fillna(method='ffill')
        .astype(int)
    )

    # Prepare data for potential TSV export
    if write_tsv_path is not None:
        export_path = (
            write_tsv_path
            if os.path.isabs(write_tsv_path)
            else os.path.join(output_dir, write_tsv_path)
        )

        (celiac_samples
            .rename_axis("Month")
            .to_frame(name="accumulative_celiac_samples")
            .assign(Month=lambda df: df.index.strftime("%Y-%m"))
            .set_index("Month")
            .to_csv(export_path, sep="\t")
        )

    # Create the plot – solid area with vertical jumps (step plot)
    plt.figure(figsize=(9, 6))
    plt.fill_between(full_months, celiac_samples.values, step='pre', alpha=0.4, color='tab:blue')

    # Set x-axis ticks to show every year
    start_year = full_months[0].year
    end_year = full_months[-1].year
    yearly_ticks = pd.date_range(start=f'{start_year}-01-01',
                                 end=f'{end_year}-12-31',
                                 freq='YS')  # YS means year start
    plt.xticks(yearly_ticks, [d.strftime('%Y') for d in yearly_ticks], rotation=45)

    # Customize the plot
    plt.title('Cumulative Number of Publically Available Celiac Samples Over Time', fontsize=12)
    plt.xlabel('Publication Date')
    plt.ylabel('Number of Available Celiac Samples')
    plt.grid(True, linestyle='--', alpha=0.7)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    plot_filename = "celiac_samples_over_time.png"
    plt.savefig(os.path.join(output_dir, plot_filename))
    plt.close()



def plot_celiac_samples_across_world_map(all_samples_df, output_dir, crop_sides=0.1, write_tsv_path=None):
    """
    Create a world map showing celiac samples per country as dots.
    The size of each dot is proportional to the number of samples.
    """
    # Add country name mapping
    country_name_map = {
        'USA': 'United States of America',
        'UK': 'United Kingdom',
        'England': 'United Kingdom',
        'U.K.': 'United Kingdom',
        'U.S.A.': 'United States of America',
        'United States': 'United States of America',
    }
    
    # Check for local map data first
    map_file = os.path.join(output_dir, "ne_110m_admin_0_countries.zip")
    
    if not os.path.exists(map_file):
        # Download world map data if not present
        import urllib.request
        print("Downloading world map data...")
        urllib.request.urlretrieve(
            "https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip",
            map_file
        )
    
    # Load world map data from local file
    world = gpd.read_file(map_file)
    
    # Count celiac samples per country with name mapping
    celiac_counts = all_samples_df[all_samples_df['Diagnosed_Celiac'] == True]['Country'].map(
        lambda x: country_name_map.get(x, x)
    ).value_counts()
    
    # Create a dictionary of country centroids
    country_centroids = {
        name: (geom.centroid.x, geom.centroid.y) 
        for name, geom in zip(world['NAME'], world['geometry'])
    }
    
    # Track missed countries
    missed_countries = []
    points_data = []
    for country, count in celiac_counts.items():
        # Find the matching country name in the world dataset
        matching_country = world[world['NAME'].str.contains(country, case=False, na=False)]
        if not matching_country.empty:
            centroid = country_centroids[matching_country.iloc[0]['NAME']]
            points_data.append({
                'Country': country,
                'Count': count,
                'geometry': Point(centroid)
            })
        else:
            missed_countries.append((country, count))
    
    # Print missed countries if any
    if missed_countries:
        print("\nWarning: Could not map the following countries:")
        for country, count in missed_countries:
            print(f"  - {country} ({count} samples)")
    
    # Create GeoDataFrame for points
    points_gdf = gpd.GeoDataFrame(points_data)
    
    # Optionally export the data used for plotting
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        (pd.DataFrame(points_data)[["Country", "Count"]]
            .to_csv(export_path, sep="\t", index=False))

    # Create the plot
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Plot world map
    world.plot(ax=ax, color='lightgrey', edgecolor='white')
    
    # Plot points
    points_gdf.plot(ax=ax, 
                   markersize=points_gdf['Count'] * 7,
                   color='red', 
                   alpha=0.4)
    
    # Customize the plot
    plt.title('Distribution of Celiac Samples Across Countries', fontsize=12)
    
    # Remove axes
    ax.axis('off')
    
    # Add legend for point sizes
    legend_sizes = [50, 100, 200, 400, 0]
    legend_elements = [plt.scatter([], [], s=size*7, 
                                 c='red', alpha=0.4, 
                                 label=f'   {size} samples ' if size > 0 else ' ')
                      for size in legend_sizes]
    ax.legend(handles=legend_elements, 
             title='Number of Samples',
             loc='lower left',
             bbox_to_anchor=(0.13, 0.175),  # Move legend away from corner
             fontsize=11,
             title_fontsize=12,
             markerscale=1.0,
             labelspacing=2.0)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "celiac_samples_world_map_temp.png"))
    plt.close()
    
    if crop_sides != False:
        if crop_sides == True:
            crop_amount = 0.1
        else:
            crop_amount = crop_sides
        # Open the saved image
        img_path = os.path.join(output_dir, "celiac_samples_world_map_temp.png")
        img = Image.open(img_path)
        
        # Calculate crop dimensions
        width, height = img.size
        crop_amount = int(width * crop_amount)
        
        # Crop the image
        cropped_img = img.crop((crop_amount, 0, width - crop_amount, height))
        
        # Save the cropped image and remove the temporary file
        final_path = os.path.join(output_dir, "celiac_samples_world_map.png")
        cropped_img.save(final_path)
        os.remove(img_path)



def plot_celiac_samples_per_sample_site(all_samples_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the number of celiac samples per sample site.
    Only includes diagnosed celiac samples, sorted in decreasing order.
    """
    # Filter for diagnosed celiac samples and count by sample site
    site_counts = all_samples_df[all_samples_df['Diagnosed_Celiac'] == True]['Sample_Site'].value_counts()

    # Capitalize the first letter of each site
    site_counts.index = site_counts.index.str.capitalize()
    
    # Optionally export the data used for plotting
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        (site_counts.to_frame(name="Count")
            .rename_axis("Sample_Site")
            .to_csv(export_path, sep="\t"))

    # Create the plot
    plt.figure(figsize=(10, 6))
    site_counts.plot(kind='bar')
    
    # Customize the plot
    plt.title('Number of Celiac Samples by Sample Site', fontsize=12)
    plt.xlabel('Sample Site')
    plt.ylabel('Number of Samples')
    plt.xticks(rotation=0)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on top of each bar
    for i, v in enumerate(site_counts):
        plt.text(i, v, str(v), ha='center', va='bottom')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "celiac_samples_per_site.png"))
    plt.close()



def plot_dataset_types(included_datasets_df, output_dir, write_tsv_path=None):
    """
    Create a table of dataset types (16S vs Shotgun, Prospective vs Non-prospective).
    Exports the results as both a TSV file and a PNG visualization.
    """
    # Create counts
    prospective = included_datasets_df['Prospective_Study'] == True
    is_16s = included_datasets_df['Sequencing_Type'] == '16S'
    
    # Initialize the table data
    table_data = {
        '': ['Non-prospective', 'Prospective'],
        '16S': [
            sum(~prospective & is_16s),
            sum(prospective & is_16s)
        ],
        'Shotgun': [
            sum(~prospective & ~is_16s),
            sum(prospective & ~is_16s)
        ]
    }
    
    # Create DataFrame
    table_df = pd.DataFrame(table_data)
    table_df = table_df.set_index('')

    # Optionally export to user-specified path
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        table_df.to_csv(export_path, sep='\t')

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Create heatmap with square cells and white-to-blue colormap
    sns.heatmap(table_df, annot=table_df.values, fmt='d', 
                cmap='Blues',  # White to blue colormap
                cbar=False, linewidths=1, linecolor='black',
                annot_kws={'size': 10, 'weight': 'normal'},
                square=True)  # Make cells square
    
    # Customize the plot
    plt.title('Dataset Types Overview', pad=10, size=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=0, ha='center')
    plt.yticks(rotation=0)
    
    # Make axis labels bold
    ax.set_xticklabels(table_df.columns)
    ax.set_yticklabels(table_df.index)
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "dataset_types_table.png"), 
                bbox_inches='tight', 
                dpi=300,
                facecolor='white')
    plt.close()



def plot_non_prospective_sample_types(all_samples_df, output_dir, write_tsv_path=None):
    """
    Create a table showing the distribution of (non-prospective) samples across diet and disease status.
    Exports the results as both a TSV file and a PNG visualization.
    Only includes samples with explicit True/False values.
    """
    # Filter for samples with valid True/False values
    valid_samples = all_samples_df[
        all_samples_df['Diagnosed_Celiac'].isin([True, False]) & 
        all_samples_df['Gluten_Free_Diet'].isin([True, False])
    ]
    
    # Create counts
    celiac = valid_samples['Diagnosed_Celiac'] == True
    gfd = valid_samples['Gluten_Free_Diet'] == True
    
    # Initialize the table data
    table_data = {
        '': ['Gluten-free', 'Normal Diet'],
        'Celiac': [
            sum(celiac & gfd),
            sum(celiac & ~gfd)
        ],
        'Healthy': [
            sum(~celiac & gfd),
            sum(~celiac & ~gfd)
        ]
    }
    
    # Create DataFrame
    table_df = pd.DataFrame(table_data)
    table_df = table_df.set_index('')

    # Optionally export to user-specified path
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        table_df.to_csv(export_path, sep='\t')

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Create heatmap with square cells and white-to-blue colormap
    sns.heatmap(table_df, annot=table_df.values, fmt='d', 
                cmap='Blues',  # White to blue colormap
                cbar=False, linewidths=1, linecolor='black',
                annot_kws={'size': 10, 'weight': 'normal'},
                square=True)  # Make cells square
    
    # Customize the plot
    plt.title('Non-Prospective Sample Types Overview', pad=10, size=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=0, ha='center')
    plt.yticks(rotation=0)
    
    # Make axis labels bold
    ax.set_xticklabels(table_df.columns)
    ax.set_yticklabels(table_df.index)
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "non_prospective_sample_types_table.png"), 
                bbox_inches='tight', 
                dpi=300,
                facecolor='white')
    plt.close()


def plot_datasets_by_amplicon_region(included_datasets_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the number of datasets per amplicon region.
    Replaces NA with "Shotgun" and sorts regions alphabetically.
    """
    # Get amplicon region counts, replacing NA with "Shotgun"
    region_counts = included_datasets_df['Amplicon_Region'].fillna('Shotgun').value_counts()
    
    # Sort alphabetically
    region_counts = region_counts.sort_index()
    
    # Optionally export the data used for plotting
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        (region_counts.to_frame(name="Count")
            .rename_axis("Amplicon_Region")
            .to_csv(export_path, sep='\t'))

    # Create the plot
    plt.figure(figsize=(10, 6))
    bars = plt.bar(region_counts.index, region_counts.values)
    
    # Customize the plot
    plt.title('Number of Datasets by Amplicon Region', fontsize=12)
    plt.xlabel('Amplicon Region')
    plt.ylabel('Number of Datasets')
    plt.xticks(rotation=0)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "datasets_by_amplicon_region.png"))
    plt.close()




def plot_prospective_sample_types(all_samples_df, output_dir, write_tsv_path=None):
    """
    Create a table showing the distribution of prospective samples based on outcome.
    Exports the results as both a TSV file and a PNG visualization.
    Only includes samples with explicit True/False values for 'Will_Develop_Celiac'.
    """
    # Filter for samples with valid True/False values for the prospective outcome
    prospective_samples = all_samples_df[all_samples_df['Will_Develop_Celiac'].isin([True, False])].copy()

    # If there are no samples, exit gracefully
    if prospective_samples.empty:
        print("Warning: No prospective samples with 'Will_Develop_Celiac' labels found. Skipping plot.")
        return

    # Create a contingency table of counts
    table_df = pd.crosstab(
        index=prospective_samples['Dataset_ID'],
        columns=prospective_samples['Will_Develop_Celiac']
    )

    # Rename columns for clarity as requested
    table_df = table_df.rename(columns={True: 'Celiac', False: 'Healthy'})

    # Ensure both 'Celiac' and 'Healthy' columns exist, filling with 0 if not
    if 'Celiac' not in table_df.columns:
        table_df['Celiac'] = 0
    if 'Healthy' not in table_df.columns:
        table_df['Healthy'] = 0
    
    # Ensure a consistent column order
    table_df = table_df[['Healthy', 'Celiac']]

    # Optionally export to user-specified path
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        table_df.to_csv(export_path, sep='\t')

    # Create figure and axis, with dynamic height
    num_datasets = len(table_df)
    fig_height = max(4, num_datasets * 0.5)
    fig, ax = plt.subplots(figsize=(8, fig_height))
    
    # Create heatmap with square cells and white-to-blue colormap
    sns.heatmap(table_df, annot=True, fmt='d', 
                cmap='Blues',
                cbar=False, linewidths=1, linecolor='black',
                annot_kws={'size': 10, 'weight': 'normal'},
                square=True)
    
    # Customize the plot
    plt.title('Prospective Sample Types', pad=10, size=12)
    plt.xticks(rotation=0, ha='center')
    plt.yticks(rotation=0)
    plt.xlabel('')
    plt.ylabel('Dataset ID')
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "prospective_sample_types_table.png"), 
                bbox_inches='tight', 
                dpi=300,
                facecolor='white')
    plt.close()


    