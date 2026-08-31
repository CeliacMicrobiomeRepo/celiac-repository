# Import libraries
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import seaborn as sns
import os
from PIL import Image
from scipy import stats

AGE_BUCKETS = {
    'Infants (0-1)': (0, 1),
    'Children (1-12)': (1, 12),
    'Adolescents (13-17)': (12, 17),
    'Adults (18-64)': (17, 64),
    'Seniors (65+)': (64, 1000)
}

def derive_raw_path(write_tsv_path, output_dir, fallback_basename):
    """
    Derive the path for a raw TSV file based on the summary TSV path.
    
    If write_tsv_path is None: return os.path.join(output_dir, fallback_basename + "_raw.tsv")
    If write_tsv_path is absolute: replace trailing ".tsv" with "_raw.tsv"
    If write_tsv_path is relative: join output_dir and replace trailing ".tsv" with "_raw.tsv"
    If it does not end in ".tsv", just append "_raw.tsv"
    """
    if write_tsv_path is None:
        return os.path.join(output_dir, fallback_basename + "_raw.tsv")
    
    if os.path.isabs(write_tsv_path):
        if write_tsv_path.endswith(".tsv"):
            return write_tsv_path[:-4] + "_raw.tsv"
        else:
            return write_tsv_path + "_raw.tsv"
    else:
        full_path = os.path.join(output_dir, write_tsv_path)
        if full_path.endswith(".tsv"):
            return full_path[:-4] + "_raw.tsv"
        else:
            return full_path + "_raw.tsv"

def write_raw_tsv(df, raw_path):
    """
    Write a dataframe to a TSV file without the index.
    """
    df.to_csv(raw_path, sep="\t", index=False)

def get_celiac_accumulation_rate(df, output_path=None):
    """
    Takes the celiac samples dataframe and calculates the growth rate
    for the 10-year period between July 2015 and June 2025.

    The results will be written to a TXT file celiac_samples_over_time.txt
    """
    
    # We first need to filter for the specific 10-year window
    # This assumes the 'Month' column is in 'YYYY-MM' string format
    start_date = '2015-07'
    end_date = '2025-06'
    
    # If Month is the index, reset it to make it a column
    if 'Month' not in df.columns and df.index.name == 'Month':
        df = df.reset_index()
    
    mask = (df['Month'] >= start_date) & (df['Month'] <= end_date)
    df_window = df.loc[mask].copy().sort_values('Month')
    
    # To run a regression, we need a numeric X axis
    # We'll just create a count of months starting from 0
    df_window['month_index'] = range(len(df_window))
    
    # The heavy lifting is done by linregress
    # It returns the slope, intercept, and statistical significance
    slope, intercept, r_val, p_val, std_err = stats.linregress(
        df_window['month_index'], 
        df_window['accumulative_celiac_samples']
    )
    
    # Pack it all into a dictionary for easy access
    results = {
        'slope_per_month': slope,
        'slope_per_year': slope * 12,
        'intercept': intercept,
        'r_squared': r_val**2,
        'p_value': p_val,
        'std_error': std_err,
        'sample_size_months': len(df_window)
    }

    if output_path:
        with open(output_path, 'w') as f:
            f.write("Celiac Accumulation Rate Analysis (July 2015 - June 2025)\n")
            f.write("========================================================\n\n")
            f.write(f"Slope (samples per month): {results['slope_per_month']:.2f}\n")
            f.write(f"Slope (samples per year):  {results['slope_per_year']:.2f}\n")
            f.write(f"Intercept:                 {results['intercept']:.2f}\n")
            f.write(f"R-squared:                {results['r_squared']:.4f}\n")
            f.write(f"P-value:                  {results['p_value']:.4e}\n")
            f.write(f"Standard Error:           {results['std_error']:.4f}\n")
            f.write(f"Sample Size (months):     {results['sample_size_months']}\n")
    
    return results


def plot_celiac_samples_over_time(all_samples_df, output_dir, write_tsv_path=None, plot_filename="celiac_samples_over_time.png", trend_results=None, include_prospective=False):
    """
    Number of celiac samples available over time.
    If include_prospective is True, includes samples where Will_Develop_Celiac is True.

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
    if include_prospective:
        mask = (all_samples_df['Diagnosed_Celiac'] == True) | (all_samples_df['Will_Develop_Celiac'] == True)
    else:
        mask = all_samples_df['Diagnosed_Celiac'] == True

    celiac_samples = (
        all_samples_df[mask]
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
        .ffill()
        .astype(int)
    )

    # Prepare data for potential TSV export and return
    celiac_samples_df = (celiac_samples
        .rename_axis("Month")
        .to_frame(name="accumulative_celiac_samples")
        .assign(
            monthly_increase_celiac_samples=lambda df: df["accumulative_celiac_samples"].diff().fillna(0).astype(int),
            Month=lambda df: df.index.strftime("%Y-%m")
        )
        .set_index("Month")
    )

    if write_tsv_path is not None:
        export_path = (
            write_tsv_path
            if os.path.isabs(write_tsv_path)
            else os.path.join(output_dir, write_tsv_path)
        )
        celiac_samples_df.to_csv(export_path, sep="\t")

    # Create the plot – solid area with vertical jumps (step plot)
    plt.figure(figsize=(9, 6))
    plt.fill_between(full_months, celiac_samples.values, step='pre', alpha=0.4, color='tab:blue')

    # If trend results are provided, plot the regression line
    if trend_results is not None:
        # We need to find the months that correspond to the regression window
        start_date = '2015-07'
        end_date = '2025-06'
        
        # Get the subset of full_months that falls within the window
        mask = (full_months >= pd.to_datetime(start_date)) & (full_months <= pd.to_datetime(end_date))
        trend_months = full_months[mask]
        
        # Calculate trend values
        # month_index starts from 0 for the first month in the window
        x = range(len(trend_months))
        y_trend = [trend_results['slope_per_month'] * i + trend_results['intercept'] for i in x]
        
        plt.plot(trend_months, y_trend, color='red', linestyle='--', label=f"Trend (R²={trend_results['r_squared']:.2f})")
        plt.legend()

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
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--', alpha=0.7)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, plot_filename))
    plt.close()

    return celiac_samples_df



def plot_celiac_samples_across_world_map(all_samples_df, output_dir, crop_sides=0.1, write_tsv_path=None, include_prospective=False):
    """
    Create a world map showing celiac samples per country as dots.
    If include_prospective is True, includes samples where Will_Develop_Celiac is True.
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
    if include_prospective:
        mask = (all_samples_df['Diagnosed_Celiac'] == True) | (all_samples_df['Will_Develop_Celiac'] == True)
    else:
        mask = all_samples_df['Diagnosed_Celiac'] == True

    celiac_counts = all_samples_df[mask]['Country'].map(
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


def plot_datasets_across_world_map(included_datasets_df, output_dir, crop_sides=0.1, write_tsv_path=None):
    """
    Create a world map showing datasets per country as dots.
    The size of each dot is proportional to the number of datasets.
    Datasets with multiple countries are counted for each country.
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
    
    # Count datasets per country, handling multiple countries per dataset
    # Split Country column by '|' and explode
    df_countries = included_datasets_df['Country'].fillna('Unknown').str.split('|')
    all_countries = [country.strip() for sublist in df_countries for country in sublist]
    
    # Map country names and count
    dataset_counts = pd.Series(all_countries).map(
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
    for country, count in dataset_counts.items():
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
            print(f"  - {country} ({count} datasets)")
    
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
    # Scale markersize for datasets (fewer datasets than samples, so maybe larger multiplier)
    points_gdf.plot(ax=ax, 
                   markersize=points_gdf['Count'] * 100,
                   color='blue', 
                   alpha=0.4)
    
    # Customize the plot
    plt.title('Distribution of Datasets Across Countries', fontsize=12)
    
    # Remove axes
    ax.axis('off')
    
    # Add legend for point sizes
    # Adjust legend sizes for datasets
    max_count = points_gdf['Count'].max()
    if max_count <= 5:
        legend_sizes = [1, 2, 5, 0]
    elif max_count <= 10:
        legend_sizes = [1, 5, 10, 0]
    else:
        legend_sizes = [1, 5, 10, 20, 0]

    legend_elements = [plt.scatter([], [], s=size*100, 
                                 c='blue', alpha=0.4, 
                                 label=f'   {size} datasets ' if size > 0 else ' ')
                      for size in legend_sizes]
    ax.legend(handles=legend_elements, 
             title='Number of Datasets',
             loc='lower left',
             bbox_to_anchor=(0.13, 0.175),
             fontsize=11,
             title_fontsize=12,
             markerscale=1.0,
             labelspacing=2.0)
    
    # Adjust layout and save
    plt.tight_layout()
    temp_path = os.path.join(output_dir, "datasets_world_map_temp.png")
    plt.savefig(temp_path)
    plt.close()
    
    if crop_sides != False:
        if crop_sides == True:
            crop_amount = 0.1
        else:
            crop_amount = crop_sides
        # Open the saved image
        img = Image.open(temp_path)
        
        # Calculate crop dimensions
        width, height = img.size
        crop_amount_px = int(width * crop_amount)
        
        # Crop the image
        cropped_img = img.crop((crop_amount_px, 0, width - crop_amount_px, height))
        
        # Save the cropped image and remove the temporary file
        final_path = os.path.join(output_dir, "datasets_world_map.png")
        cropped_img.save(final_path)
        os.remove(temp_path)


def plot_celiac_samples_per_sample_site(all_samples_df, output_dir, write_tsv_path=None, include_prospective=False):
    """
    Create a bar plot showing the number of celiac samples per sample site.
    If include_prospective is True, includes samples where Will_Develop_Celiac is True.
    Sorted in decreasing order.
    """
    # Filter for celiac samples and count by sample site
    if include_prospective:
        mask = (all_samples_df['Diagnosed_Celiac'] == True) | (all_samples_df['Will_Develop_Celiac'] == True)
    else:
        mask = all_samples_df['Diagnosed_Celiac'] == True
        
    site_counts = all_samples_df[mask]['Sample_Site'].value_counts()

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


def plot_datasets_by_site(included_datasets_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the number of datasets per sample site.
    Datasets with multiple sites are counted for each site.
    Sorted in decreasing order.
    """
    # Split Sample_Sites and explode to count each site individually
    # First, handle potential NA values
    df_sites = included_datasets_df['Sample_Sites'].fillna('Unknown').str.split('|')
    
    # Flatten the list of sites and count them
    all_sites = [site for sublist in df_sites for site in sublist]
    site_counts = pd.Series(all_sites).value_counts()

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
    plt.title('Number of Datasets by Sample Site', fontsize=12)
    plt.xlabel('Sample Site')
    plt.ylabel('Number of Datasets')
    plt.xticks(rotation=0)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on top of each bar
    for i, v in enumerate(site_counts):
        plt.text(i, v, str(v), ha='center', va='bottom')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "datasets_by_site.png"))
    plt.close()


def plot_samples_per_dataset_significant_factors(all_samples_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the number of samples per dataset with, without, or
    with unresolved significant factors.
    Exports the results as both a TSV file and a PNG visualization.
    """
    df_sig = all_samples_df.copy()
    df_sig['Any_Significant_Factor'] = (
        df_sig['Any_Significant_Factor']
        .map({True: 'With', False: 'Without', 'True': 'With', 'False': 'Without'})
        .fillna('Unknown')
    )
    
    # Create a contingency table of counts
    table_df = pd.crosstab(
        index=df_sig['Dataset_ID'],
        columns=df_sig['Any_Significant_Factor']
    )

    # Ensure all three states are represented, even when a state has no samples.
    for factor_state in ['Without', 'With', 'Unknown']:
        if factor_state not in table_df.columns:
            table_df[factor_state] = 0
    
    # Without is at the bottom of the stacked bar; unresolved values are on top.
    table_df = table_df[['Without', 'With', 'Unknown']]

    # Prepare DataFrame for TSV export with requested column names
    export_df = table_df.copy().reset_index()
    export_df.columns = [
        'Dataset_ID',
        'Num_Samples_Without_Significant_Factor',
        'Num_Samples_With_Significant_Factor',
        'Num_Samples_With_Unknown_Significant_Factor_Status'
    ]
    
    export_df = export_df[[
        'Dataset_ID',
        'Num_Samples_With_Significant_Factor',
        'Num_Samples_Without_Significant_Factor',
        'Num_Samples_With_Unknown_Significant_Factor_Status'
    ]]

    # Optionally export to user-specified path
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        export_df.to_csv(export_path, sep='\t', index=False)

    # Create the plot
    plt.figure(figsize=(max(10, len(table_df) * 0.5), 6))
    
    # Plot stacked bars: blue for Without, red for With, and grey for Unknown.
    table_df.plot(
        kind='bar',
        stacked=True,
        color=['blue', 'red', 'grey'],
        width=0.8,
        figsize=(max(10, len(table_df) * 0.5), 6)
    )
    
    # Customize the plot
    plt.title('Number of Samples per Dataset by Significant Factors', fontsize=12)
    plt.xlabel('Dataset ID')
    plt.ylabel('Number of Samples')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.legend(
        ['No significant factors', 'With significant factors', 'Unknown / unresolved'],
        title='Significant Factors'
    )
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "samples_per_dataset_significant_factors.png"))
    plt.close()


def plot_samples_per_group(all_samples_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the number of samples per group.
    Groups are: HC, ACD, HC_GFD, TCD, PHC, PCD.
    Full names are used for the x-axis labels.
    Exports the results as both a TSV file and a PNG visualization.
    """
    # Define the group mapping and order
    group_mapping = {
        'HC': 'Healthy Control',
        'ACD': 'Active Celiac',
        'HC_GFD': 'Healthy Control on a GFD',
        'TCD': 'Treated Celiac',
        'PHC': 'Prospective Healthy Control',
        'PCD': 'Prospective Celiac Disease'
    }
    group_order = ['HC', 'ACD', 'HC_GFD', 'TCD', 'PHC', 'PCD']
    
    # Count the number of samples per group
    group_counts = all_samples_df['Group'].value_counts()
    
    # Create a DataFrame for the counts in the specified order
    plot_df = pd.DataFrame({
        'Group': group_order,
        'Num_Samples': [group_counts.get(group, 0) for group in group_order],
        'Full_Name': [group_mapping[group] for group in group_order]
    })
    
    # Optionally export the data used for plotting
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        plot_df[['Group', 'Num_Samples']].to_csv(export_path, sep='\t', index=False)

    # Create the plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(plot_df['Full_Name'], plot_df['Num_Samples'], color='skyblue', edgecolor='black')
    
    # Customize the plot
    plt.title('Number of Samples per Group', fontsize=14)
    plt.xlabel('Group', fontsize=12)
    plt.ylabel('Number of Samples', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 5,
                 f'{int(height)}',
                 ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "samples_per_group.png"), dpi=300)
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


def analyse_sex(all_samples_df, included_datasets_df, output_dir, write_tsv_path="sex_proportions.tsv"):
    """
    Analyse the distribution of sex for datasets with sample-level sex metadata.
    """
    datasets_with_sex = included_datasets_df[
        included_datasets_df['Sex_Metadata_Availability'] == 'sample_level'
    ]['Dataset_ID']
    
    if datasets_with_sex.empty:
        print("No datasets with sex metadata found.")
        return

    # Filter all samples for these datasets
    df_sex = all_samples_df[all_samples_df['Dataset_ID'].isin(datasets_with_sex)].copy()
    
    # Normalize sex column to lower case
    df_sex['Sex'] = df_sex['Sex'].astype(str).str.lower()
    
    # Calculate proportions per dataset
    results = []
    for dataset_id in datasets_with_sex:
        subset = df_sex[df_sex['Dataset_ID'] == dataset_id]
        total = len(subset)
        if total == 0:
            continue
            
        male_count = len(subset[subset['Sex'] == 'male'])
        female_count = len(subset[subset['Sex'] == 'female'])
        
        # Calculate proportions based on total samples, assuming remaining are unknown/other
        # However, prompt asks for proportion of samples that are male and female.
        # Usually implies denominator is total samples.
        prop_male = male_count / total
        prop_female = female_count / total
        
        results.append({
            'Dataset_ID': dataset_id,
            'Sex_Proportion_Female': prop_female,
            'Sex_Proportion_Male': prop_male
        })
    
    results_df = pd.DataFrame(results)
    
    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        results_df.to_csv(export_path, sep='\t', index=False)
        
    # Record stats
    stats_path = os.path.join(output_dir, "sex_proportions.txt")
    with open(stats_path, 'w') as f:
        f.write("Sex Proportion Analysis\n")
        f.write("=======================\n\n")
        
        for col in ['Sex_Proportion_Female', 'Sex_Proportion_Male']:
            f.write(f"{col}:\n")
            f.write(f"  Min:  {results_df[col].min():.4f}\n")
            f.write(f"  Mean: {results_df[col].mean():.4f}\n")
            f.write(f"  Max:  {results_df[col].max():.4f}\n\n")

        # Pooled sample-level statistic across ALL datasets (not just sample_level ones)
        sex_all = all_samples_df['Sex'].astype(str).str.lower()
        n_female = (sex_all == 'female').sum()
        n_male = (sex_all == 'male').sum()
        n_known = n_female + n_male
        f.write("Pooled (all samples with known sex):\n")
        f.write(f"  Female: {n_female}/{n_known} = {100 * n_female / n_known:.1f}%\n")
        f.write(f"  Male:   {n_male}/{n_known} = {100 * n_male / n_known:.1f}%\n\n")

    # Plot stacked bar chart
    if not results_df.empty:
        plot_df = results_df.set_index('Dataset_ID')
        
        plt.figure(figsize=(12, 6))
        
        # Plot bars
        # Blue for male, Red for female
        # We can use plot(kind='bar', stacked=True)
        # We need a dataframe with columns 'Male' and 'Female' for plotting
        plot_data = plot_df[['Sex_Proportion_Male', 'Sex_Proportion_Female']].rename(
            columns={'Sex_Proportion_Male': 'Male', 'Sex_Proportion_Female': 'Female'}
        )
        
        ax = plot_data.plot(kind='bar', stacked=True, color=['blue', 'red'], width=0.8, figsize=(12, 6))
        
        plt.title('Sex Proportions per Dataset', fontsize=12)
        plt.xlabel('Dataset ID')
        plt.ylabel('Sex Proportion')
        plt.ylim(0, 1.0)
        plt.grid(True, axis='y', linestyle='--', alpha=0.7)
        plt.legend(title='Sex')
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "sex_proportions.png"))
        plt.close()


def analyse_age(all_samples_df, included_datasets_df, output_dir, write_tsv_path="age_distribution.tsv"):
    """
    Analyse the distribution of age for all datasets.
    """
    results = []
    
    for _, row in included_datasets_df.iterrows():
        dataset_id = row['Dataset_ID']
        has_sample_level_age = row['Age_Metadata_Availability'] == 'sample_level'
        age_range_str = str(row['Age_Range'])
        
        age_min = pd.NA
        age_mean = pd.NA
        age_max = pd.NA
        age_sd = pd.NA
        
        # Parse Age_Range for min/max fallback
        range_min = pd.NA
        range_max = pd.NA
        
        try:
            if pd.notna(age_range_str) and age_range_str.lower() != 'nan' and '-' in age_range_str:
                parts = age_range_str.split('-')
                if len(parts) == 2:
                    range_min = float(parts[0])
                    range_max = float(parts[1])
        except ValueError:
            pass # Keep as NA if parsing fails
            
        
        if has_sample_level_age:
            # Get samples for this dataset
            samples = all_samples_df[all_samples_df['Dataset_ID'] == dataset_id]
            # Convert Age to numeric, coercing errors to NaN
            ages = pd.to_numeric(samples['Age'], errors='coerce').dropna()
            
            if not ages.empty:
                age_min = ages.min()
                age_max = ages.max()
                age_mean = ages.mean()
                age_sd = ages.std()
                
                # Check for discrepancies with Age_Range
                if pd.notna(range_min) and pd.notna(range_max):
                    # Allow for small floating point differences
                    if (age_min < range_min - 0.1) or (age_max > range_max + 0.1):
                        print(f"Warning: Age discrepancy in {dataset_id}. Samples range: {age_min}-{age_max}, Metadata range: {range_min}-{range_max}")
            else:
                # If no valid ages found in samples, fallback to range
                age_min = range_min
                age_max = range_max
                age_mean = "NA"
                age_sd = "NA"
                
        else:
            # Use Age_Range
            age_min = range_min
            age_max = range_max
            age_mean = "NA"
            age_sd = "NA"
            
        results.append({
            'Dataset_ID': dataset_id,
            'Age_Min': age_min,
            'Age_Mean': age_mean,
            'Age_Max': age_max,
            'Age_SD': age_sd
        })
        
    results_df = pd.DataFrame(results)
    
    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        results_df.to_csv(export_path, sep='\t', index=False)
        
    # Plot age ranges
    # We need Min and Max to be numeric for plotting
    plot_df = results_df.dropna(subset=['Age_Min', 'Age_Max']).copy()
    plot_df['Age_Min'] = pd.to_numeric(plot_df['Age_Min'])
    plot_df['Age_Max'] = pd.to_numeric(plot_df['Age_Max'])
    
    # Sort for better visualization (e.g. by Age_Min)
    plot_df = plot_df.sort_values('Age_Min', ascending=False)
    
    plt.figure(figsize=(10, max(6, len(plot_df) * 0.4)))
    
    # Create horizontal lines
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        dataset = row['Dataset_ID']
        min_val = row['Age_Min']
        max_val = row['Age_Max']
        
        plt.hlines(y=i, xmin=min_val, xmax=max_val, linewidth=3, color='skyblue')
        # Add dots at ends
        plt.plot(min_val, i, 'o', color='skyblue', markersize=6)
        plt.plot(max_val, i, 'o', color='skyblue', markersize=6)
        
    plt.yticks(range(len(plot_df)), plot_df['Dataset_ID'])
    plt.xlabel('Age (Years)')
    plt.ylabel('Dataset ID')
    plt.title('Age Ranges of Datasets')
    plt.xlim(0, 110)
    plt.grid(True, axis='x', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "age_distribution.png"))
    plt.close()


def analyse_read_retention_16S(all_samples_df, output_dir, write_tsv_path="read_retention_16S.tsv"):
    """
    Analyse read retention over steps for 16S datasets.
    """
    # Filter for 16S datasets
    df_16s = all_samples_df[all_samples_df['Sequencing_Type'] == '16S'].copy()

    if df_16s.empty:
        print("No 16S datasets found for read retention analysis.")
        return

    # Calculate percentages
    # Ensure numeric types
    cols = ['Num_Reads_Input', 'Num_Reads_Filtered', 'Num_Reads_DenoisedF', 'Num_Reads_Nonchim']
    for col in cols:
        df_16s[col] = pd.to_numeric(df_16s[col], errors='coerce')

    # Drop rows with NaN in Input if any, to avoid division by zero or NaN results
    df_16s = df_16s.dropna(subset=['Num_Reads_Input'])
    # Avoid division by zero
    df_16s = df_16s[df_16s['Num_Reads_Input'] > 0]

    df_16s['Filtered_Percent_Reads_Remaining'] = (df_16s['Num_Reads_Filtered'] / df_16s['Num_Reads_Input']) * 100
    df_16s['DenoisedF_Percent_Reads_Remaining'] = (df_16s['Num_Reads_DenoisedF'] / df_16s['Num_Reads_Input']) * 100
    df_16s['Nonchim_Percent_Reads_Remaining'] = (df_16s['Num_Reads_Nonchim'] / df_16s['Num_Reads_Input']) * 100

    # Group by Dataset_ID and calculate stats
    stats_cols = [
        'Filtered_Percent_Reads_Remaining', 
        'DenoisedF_Percent_Reads_Remaining', 
        'Nonchim_Percent_Reads_Remaining'
    ]
    
    agg_funcs = ['min', 'mean', 'median', 'max', 'std']
    
    grouped = df_16s.groupby('Dataset_ID')[stats_cols].agg(agg_funcs)
    
    # Flatten MultiIndex columns and rename
    new_columns = []
    for col_name, stat_name in grouped.columns:
        if stat_name == 'std':
            stat_suffix = 'SD'
        else:
            stat_suffix = stat_name.capitalize()
        new_columns.append(f"{col_name}_{stat_suffix}")
        
    grouped.columns = new_columns
    grouped = grouped.reset_index()
    
    # Reorder columns to match requirement
    ordered_cols = ['Dataset_ID']
    for base_col in stats_cols:
        for stat in ['Min', 'Mean', 'Median', 'Max', 'SD']:
            col_name = f"{base_col}_{stat}"
            if col_name in grouped.columns:
                ordered_cols.append(col_name)
            
    final_stats_df = grouped[ordered_cols]
    
    # Add an "Overall" row with aggregate stats across datasets
    overall_row = {'Dataset_ID': 'Overall'}
    for col in final_stats_df.columns:
        if col == 'Dataset_ID':
            continue
        if col.endswith('_Min'):
            overall_row[col] = final_stats_df[col].min()
        elif col.endswith('_Max'):
            overall_row[col] = final_stats_df[col].max()
        elif col.endswith('_Mean') or col.endswith('_SD'):
            overall_row[col] = final_stats_df[col].mean()
        elif col.endswith('_Median'):
            overall_row[col] = final_stats_df[col].median()
            
    final_stats_df = pd.concat([final_stats_df, pd.DataFrame([overall_row])], ignore_index=True)
    
    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        final_stats_df.to_csv(export_path, sep='\t', index=False)

    # Plot
    # Prepare data for plotting (melt)
    plot_df = df_16s[['Dataset_ID'] + stats_cols].melt(
        id_vars=['Dataset_ID'], 
        value_vars=stats_cols, 
        var_name='Step', 
        value_name='Percentage'
    )
    
    # Clean up Step names for legend
    step_map = {
        'Filtered_Percent_Reads_Remaining': 'Filtered',
        'DenoisedF_Percent_Reads_Remaining': 'DenoisedF',
        'Nonchim_Percent_Reads_Remaining': 'Nonchim'
    }
    plot_df['Step'] = plot_df['Step'].map(step_map)
    
    # Write raw TSV with exactly what goes into seaborn
    raw_path = derive_raw_path(write_tsv_path, output_dir, "read_retention_16S")
    write_raw_tsv(plot_df, raw_path)
    
    # Plot using seaborn boxplot
    num_datasets = plot_df['Dataset_ID'].nunique()
    plt.figure(figsize=(max(10, num_datasets * 0.5), 6))
    
    datasets = sorted(plot_df['Dataset_ID'].unique())
    
    palette = {'Filtered': 'blue', 'DenoisedF': 'green', 'Nonchim': 'red'}
    
    sns.boxplot(
        data=plot_df, 
        x='Dataset_ID', 
        y='Percentage', 
        hue='Step', 
        palette=palette,
        hue_order=['Filtered', 'DenoisedF', 'Nonchim'],
        order=datasets
    )
    
    plt.title('Read Retention over Steps for 16S Datasets')
    plt.ylabel('Percentage of Reads Remaining')
    plt.xlabel('Dataset ID')
    plt.ylim(0, 105)
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.legend(title='DADA2 Step')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "read_retention_16S.png"))
    plt.close()


def split_datasets_by_site(df):
    """
    Handle datasets with multiple sample sites by splitting them into separate entries.
    For example, if SG_80_Mouzan has stool and duodenal samples, it becomes
    SG_80_Mouzan_Stool and SG_80_Mouzan_Duodenal.
    """
    # First, calculate unique sample sites per dataset
    site_counts = df.groupby('Dataset_ID')['Sample_Site'].nunique()
    datasets_with_multiple_sites = site_counts[site_counts > 1].index
    
    if not datasets_with_multiple_sites.empty:
        # Function to generate new ID
        def get_display_id(row):
            if row['Dataset_ID'] in datasets_with_multiple_sites:
                # Clean up sample site (capitalize)
                site = str(row['Sample_Site']).capitalize() if pd.notna(row['Sample_Site']) else "Unknown"
                return f"{row['Dataset_ID']}_{site}"
            return row['Dataset_ID']
            
        df['Dataset_ID'] = df.apply(get_display_id, axis=1)
    
    return df


def analyse_host_read_removal_shotgun(all_samples_df, output_dir, write_tsv_path="host_read_removal.tsv"):
    """
    Analyse Host Read Removal for Shotgun datasets.
    """
    # Filter for Shotgun datasets (Sequencing_Type == 'SG')
    df_sg = all_samples_df[all_samples_df['Sequencing_Type'] == 'SG'].copy()

    if df_sg.empty:
        print("No Shotgun datasets found for host read removal analysis.")
        return

    # Ensure numeric
    df_sg['Percent_Host_Reads_Removed'] = pd.to_numeric(df_sg['Percent_Host_Reads_Removed'], errors='coerce')
    
    # Drop rows where Percent_Host_Reads_Removed is NaN
    df_sg = df_sg.dropna(subset=['Percent_Host_Reads_Removed'])
    
    if df_sg.empty:
        print("No Shotgun samples with valid Percent_Host_Reads_Removed found.")
        return

    # Handle datasets with multiple sample sites by splitting them
    df_sg = split_datasets_by_site(df_sg)

    # Calculate stats per dataset
    agg_funcs = ['min', 'mean', 'median', 'max', 'std']
    grouped = df_sg.groupby('Dataset_ID')['Percent_Host_Reads_Removed'].agg(agg_funcs)
    
    # Rename columns
    grouped.columns = [f"Percent_Host_Reads_Removed_{col.capitalize() if col != 'std' else 'SD'}" for col in grouped.columns]
    grouped = grouped.reset_index()
    
    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        grouped.to_csv(export_path, sep='\t', index=False)

    # Plot
    # Calculate Reads Remaining
    df_sg['Percent_Reads_Remaining'] = 100 - df_sg['Percent_Host_Reads_Removed']
    
    # Sort datasets
    datasets = sorted(df_sg['Dataset_ID'].unique())
    
    # Create plot dataframe with exactly what goes into seaborn
    plot_df = df_sg[['Dataset_ID', 'Percent_Reads_Remaining']].copy()
    
    # Write raw TSV
    raw_path = derive_raw_path(write_tsv_path, output_dir, "host_read_removal")
    write_raw_tsv(plot_df, raw_path)
    
    plt.figure(figsize=(max(10, len(datasets) * 0.5), 6))
    
    sns.boxplot(
        data=plot_df,
        x='Dataset_ID',
        y='Percent_Reads_Remaining',
        color='blue', # User asked for blue box.
        order=datasets
    )
    
    plt.title('Host Read Removal for Shotgun Datasets')
    plt.ylabel('Percentage of Reads Remaining')
    plt.xlabel('Dataset ID')
    plt.ylim(0, 100)
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', edgecolor='black', label='Percent Reads Remaining')]
    plt.legend(handles=legend_elements, loc='best', title='Measure')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "host_read_removal_shotgun.png"))
    plt.close()


def analyse_reads_remaining_after_dada2_16S(all_samples_df, output_dir, write_tsv_path="reads_remaining_after_dada2_16S.tsv"):
    """
    Analyse and plot numbers of reads remaining after DADA2 for 16S.
    """
    # Filter for 16S datasets
    df_16s = all_samples_df[all_samples_df['Sequencing_Type'] == '16S'].copy()

    if df_16s.empty:
        print("No 16S datasets found for reads remaining analysis.")
        return

    # Ensure numeric
    df_16s['Num_Reads_Nonchim'] = pd.to_numeric(df_16s['Num_Reads_Nonchim'], errors='coerce')
    
    # Drop rows where Num_Reads_Nonchim is NaN
    df_16s = df_16s.dropna(subset=['Num_Reads_Nonchim'])

    # Calculate stats per dataset
    agg_funcs = ['min', 'mean', 'median', 'max', 'std']
    grouped = df_16s.groupby('Dataset_ID')['Num_Reads_Nonchim'].agg(agg_funcs)
    
    # Rename columns
    new_columns = []
    for col in grouped.columns:
        if col == 'std':
            suffix = 'SD'
        else:
            suffix = col.capitalize()
        new_columns.append(f"Num_Reads_Nonchim_{suffix}")
    
    grouped.columns = new_columns
    grouped = grouped.reset_index()
    
    # Reorder columns to ensure correct order
    ordered_cols = ['Dataset_ID']
    for stat in ['Min', 'Mean', 'Median', 'Max', 'SD']:
        col_name = f"Num_Reads_Nonchim_{stat}"
        if col_name in grouped.columns:
            ordered_cols.append(col_name)
            
    final_stats_df = grouped[ordered_cols]
    
    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        final_stats_df.to_csv(export_path, sep='\t', index=False)

    # Plot
    # Sort datasets for consistent plotting
    datasets = sorted(df_16s['Dataset_ID'].unique())
    
    # Create plot dataframe with exactly what goes into seaborn
    plot_df = df_16s[['Dataset_ID', 'Num_Reads_Nonchim']].copy()
    
    # Write raw TSV
    raw_path = derive_raw_path(write_tsv_path, output_dir, "reads_remaining_after_dada2_16S")
    write_raw_tsv(plot_df, raw_path)
    
    plt.figure(figsize=(max(10, len(datasets) * 0.5), 6))
    
    sns.boxplot(
        data=plot_df,
        x='Dataset_ID',
        y='Num_Reads_Nonchim',
        color='blue', # User asked for blue box
        order=datasets
    )
    
    plt.title('Number of Reads Remaining After DADA2 for 16S Datasets')
    plt.ylabel('Number of Reads (Nonchim)')
    plt.xlabel('Dataset ID')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', edgecolor='black', label='Num Reads Nonchim')]
    plt.legend(handles=legend_elements, loc='best', title='Measure')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "reads_remaining_after_dada2_16S.png"))
    plt.close()


def analyse_num_asvs_16S(all_samples_df, output_dir, write_tsv_path="num_asvs.tsv"):
    """
    Analyse and plot numbers of unique ASVs per sample for 16S.
    """
    # Filter for 16S datasets
    df_16s = all_samples_df[all_samples_df['Sequencing_Type'] == '16S'].copy()

    if df_16s.empty:
        print("No 16S datasets found for ASV analysis.")
        return

    # Ensure numeric
    df_16s['Num_ASVs'] = pd.to_numeric(df_16s['Num_ASVs'], errors='coerce')
    
    # Drop rows where Num_ASVs is NaN
    df_16s = df_16s.dropna(subset=['Num_ASVs'])

    if df_16s.empty:
        print("No 16S samples with valid Num_ASVs found.")
        return

    # Calculate stats per dataset
    agg_funcs = ['min', 'mean', 'median', 'max', 'std']
    grouped = df_16s.groupby('Dataset_ID')['Num_ASVs'].agg(agg_funcs)
    
    # Rename columns
    new_columns = []
    for col in grouped.columns:
        if col == 'std':
            suffix = 'SD'
        else:
            suffix = col.capitalize()
        new_columns.append(f"Num_ASVs_{suffix}")
    
    grouped.columns = new_columns
    grouped = grouped.reset_index()
    
    # Reorder columns to ensure correct order
    ordered_cols = ['Dataset_ID']
    for stat in ['Min', 'Mean', 'Median', 'Max', 'SD']:
        col_name = f"Num_ASVs_{stat}"
        if col_name in grouped.columns:
            ordered_cols.append(col_name)
            
    final_stats_df = grouped[ordered_cols]
    
    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        final_stats_df.to_csv(export_path, sep='\t', index=False)

    # Plot
    # Sort datasets for consistent plotting
    datasets = sorted(df_16s['Dataset_ID'].unique())
    
    # Create plot dataframe with exactly what goes into seaborn
    plot_df = df_16s[['Dataset_ID', 'Num_ASVs']].copy()
    
    # Write raw TSV
    raw_path = derive_raw_path(write_tsv_path, output_dir, "num_asvs")
    write_raw_tsv(plot_df, raw_path)
    
    plt.figure(figsize=(max(10, len(datasets) * 0.5), 6))
    
    sns.boxplot(
        data=plot_df,
        x='Dataset_ID',
        y='Num_ASVs',
        color='blue', # Consistent with other plots
        order=datasets
    )
    
    plt.title('Number of Unique ASVs per Sample for 16S Datasets')
    plt.ylabel('Number of Unique ASVs')
    plt.xlabel('Dataset ID')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', edgecolor='black', label='Num ASVs')]
    plt.legend(handles=legend_elements, loc='best', title='Measure')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "num_asvs.png"))
    plt.close()


def analyse_num_sgbs_shotgun(all_samples_df, output_dir, write_tsv_path="num_sgbs.tsv"):
    """
    Analyse and plot numbers of unique SGBs per sample for Shotgun.
    """
    # Filter for Shotgun datasets (Sequencing_Type == 'SG')
    df_sg = all_samples_df[all_samples_df['Sequencing_Type'] == 'SG'].copy()

    if df_sg.empty:
        print("No Shotgun datasets found for SGB analysis.")
        return

    # Ensure numeric
    df_sg['Num_SGBs'] = pd.to_numeric(df_sg['Num_SGBs'], errors='coerce')
    
    # Drop rows where Num_SGBs is NaN
    df_sg = df_sg.dropna(subset=['Num_SGBs'])

    if df_sg.empty:
        print("No Shotgun samples with valid Num_SGBs found.")
        return

    # Handle datasets with multiple sample sites by splitting them
    df_sg = split_datasets_by_site(df_sg)

    # Calculate stats per dataset
    agg_funcs = ['min', 'mean', 'median', 'max', 'std']
    grouped = df_sg.groupby('Dataset_ID')['Num_SGBs'].agg(agg_funcs)
    
    # Rename columns
    new_columns = []
    for col in grouped.columns:
        if col == 'std':
            suffix = 'SD'
        else:
            suffix = col.capitalize()
        new_columns.append(f"Num_SGBs_{suffix}")
    
    grouped.columns = new_columns
    grouped = grouped.reset_index()
    
    # Reorder columns to ensure correct order
    ordered_cols = ['Dataset_ID']
    for stat in ['Min', 'Mean', 'Median', 'Max', 'SD']:
        col_name = f"Num_SGBs_{stat}"
        if col_name in grouped.columns:
            ordered_cols.append(col_name)
            
    final_stats_df = grouped[ordered_cols]
    
    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        final_stats_df.to_csv(export_path, sep='\t', index=False)

    # Plot
    # Sort datasets for consistent plotting
    datasets = sorted(df_sg['Dataset_ID'].unique())
    
    # Create plot dataframe with exactly what goes into seaborn
    plot_df = df_sg[['Dataset_ID', 'Num_SGBs']].copy()
    
    # Write raw TSV
    raw_path = derive_raw_path(write_tsv_path, output_dir, "num_sgbs")
    write_raw_tsv(plot_df, raw_path)
    
    plt.figure(figsize=(max(10, len(datasets) * 0.5), 6))
    
    sns.boxplot(
        data=plot_df,
        x='Dataset_ID',
        y='Num_SGBs',
        color='blue', # Consistent with other plots
        order=datasets
    )
    
    plt.title('Number of Unique SGBs per Sample for Shotgun Datasets')
    plt.ylabel('Number of Unique SGBs')
    plt.xlabel('Dataset ID')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', edgecolor='black', label='Num SGBs')]
    plt.legend(handles=legend_elements, loc='best', title='Measure')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "num_sgbs.png"))
    plt.close()


def analyse_asv_classification_16S(included_datasets_df, output_dir, write_tsv_path_absolute=None, write_tsv_path_relative=None):
    """
    Analyse and plot numbers of ASVs classified at genus and species level for 16S.
    Produces two plots:
    1. Absolute number of ASVs classified (asv_classification_16S_absolute.png)
    2. Percentage of ASVs classified relative to Total_Num_ASVs (asv_classification_16S_relative.png)
    """
    # Filter for 16S datasets
    df_16s = included_datasets_df[included_datasets_df['Sequencing_Type'] == '16S'].copy()
    
    if df_16s.empty:
        print("No 16S datasets found for ASV classification analysis.")
        return
        
    # Ensure relevant columns are numeric
    cols = ['Total_Num_ASVs', 'Num_ASVs_Classified_Genus', 'Num_ASVs_Classified_Species']
    for col in cols:
        df_16s[col] = pd.to_numeric(df_16s[col], errors='coerce')
        
    # Drop rows with NaN in these columns if any
    df_16s = df_16s.dropna(subset=cols)
    
    if df_16s.empty:
        print("No 16S datasets with valid ASV classification counts found.")
        return
        
    # Sort by Dataset_ID for consistent plotting
    df_16s = df_16s.sort_values('Dataset_ID')
    datasets = df_16s['Dataset_ID'].tolist()
    
    # --- Plot 1: Absolute Numbers ---

    # Optionally export absolute data
    if write_tsv_path_absolute is not None:
        export_path = write_tsv_path_absolute if os.path.isabs(write_tsv_path_absolute) else os.path.join(output_dir, write_tsv_path_absolute)
        df_16s[['Dataset_ID', 'Num_ASVs_Classified_Genus', 'Num_ASVs_Classified_Species']].to_csv(export_path, sep='\t', index=False)
    
    plt.figure(figsize=(max(10, len(datasets) * 0.5), 6))
    
    x = range(len(datasets))
    width = 0.35
    
    plt.bar([i - width/2 for i in x], df_16s['Num_ASVs_Classified_Genus'], width, label='Genus', color='blue')
    plt.bar([i + width/2 for i in x], df_16s['Num_ASVs_Classified_Species'], width, label='Species', color='red')
    
    plt.title('Number of ASVs Classified at Genus and Species Level (16S)')
    plt.ylabel('Number of ASVs')
    plt.xlabel('Dataset ID')
    plt.xticks(x, datasets, rotation=45, ha='right')
    plt.legend()
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "asv_classification_16S_absolute.png"))
    plt.close()
    
    # --- Plot 2: Relative Percentages ---
    
    # Calculate percentages
    df_16s['Percent_Genus'] = (df_16s['Num_ASVs_Classified_Genus'] / df_16s['Total_Num_ASVs']) * 100
    df_16s['Percent_Species'] = (df_16s['Num_ASVs_Classified_Species'] / df_16s['Total_Num_ASVs']) * 100

    # Optionally export relative data
    if write_tsv_path_relative is not None:
        export_path = write_tsv_path_relative if os.path.isabs(write_tsv_path_relative) else os.path.join(output_dir, write_tsv_path_relative)
        df_16s[['Dataset_ID', 'Percent_Genus', 'Percent_Species']].to_csv(export_path, sep='\t', index=False)
    
    plt.figure(figsize=(max(10, len(datasets) * 0.5), 6))
    
    plt.bar([i - width/2 for i in x], df_16s['Percent_Genus'], width, label='Genus', color='blue')
    plt.bar([i + width/2 for i in x], df_16s['Percent_Species'], width, label='Species', color='red')
    
    plt.title('Percentage of ASVs Classified at Genus and Species Level (16S)')
    plt.ylabel('Percentage of ASVs (%)')
    plt.xlabel('Dataset ID')
    plt.xticks(x, datasets, rotation=45, ha='right')
    plt.legend()
    plt.ylim(0, 100)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "asv_classification_16S_relative.png"))
    plt.close()


def analyse_extraction_kit_type(all_samples_df, output_dir, write_tsv_path="reads_extraction_kit_type.tsv"):
    """
    Analyse Extraction Kit type against read counts for 16S (DNA_Extraction_Is_Mechanical and Num_Reads_Nonchim).
    """
    # Filter for 16S datasets
    df_16s = all_samples_df[all_samples_df['Sequencing_Type'] == '16S'].copy()

    if df_16s.empty:
        print("No 16S datasets found for extraction kit analysis.")
        return

    # Ensure numeric
    df_16s['Num_Reads_Nonchim'] = pd.to_numeric(df_16s['Num_Reads_Nonchim'], errors='coerce')
    
    # Drop rows where Num_Reads_Nonchim is NaN
    df_16s = df_16s.dropna(subset=['Num_Reads_Nonchim'])

    # Group by Dataset_ID to get median reads and extraction type
    # We assume DNA_Extraction_Is_Mechanical is consistent within a dataset
    dataset_stats = df_16s.groupby('Dataset_ID').agg({
        'Num_Reads_Nonchim': 'median',
        'DNA_Extraction_Is_Mechanical': 'first'
    }).reset_index()
    
    dataset_stats.rename(columns={'Num_Reads_Nonchim': 'Num_Reads_Nonchim_Median'}, inplace=True)

    # Write TSV
    if write_tsv_path:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        dataset_stats.to_csv(export_path, sep='\t', index=False)

    # Prepare for plotting and stats
    # Ensure DNA_Extraction_Is_Mechanical is boolean
    if dataset_stats['DNA_Extraction_Is_Mechanical'].dtype == object:
        dataset_stats['DNA_Extraction_Is_Mechanical'] = dataset_stats['DNA_Extraction_Is_Mechanical'].map({'True': True, 'False': False, 'TRUE': True, 'FALSE': False})
    
    dataset_stats['DNA_Extraction_Is_Mechanical'] = dataset_stats['DNA_Extraction_Is_Mechanical'].astype(bool)
    
    mech_true = dataset_stats[dataset_stats['DNA_Extraction_Is_Mechanical'] == True]['Num_Reads_Nonchim_Median']
    mech_false = dataset_stats[dataset_stats['DNA_Extraction_Is_Mechanical'] == False]['Num_Reads_Nonchim_Median']

    # Statistical Test (Mann-Whitney U)
    # We use Mann-Whitney U because sample sizes might be small and distribution non-normal
    stat, p_val = stats.mannwhitneyu(mech_true, mech_false, alternative='two-sided')
    
    # Write stats to text file
    stats_path = os.path.join(output_dir, "reads_extraction_kit_type.txt")
    with open(stats_path, 'w') as f:
        f.write("Extraction Kit Type Analysis (Mechanical Lysis vs No Mechanical Lysis)\n")
        f.write("==================================================================\n\n")
        f.write(f"Mechanical Lysis (True):  n={len(mech_true)}\n")
        f.write(f"  Median: {mech_true.median():.2f}\n")
        f.write(f"  Mean:   {mech_true.mean():.2f}\n")
        f.write(f"  Std:    {mech_true.std():.2f}\n\n")
        
        f.write(f"No Mechanical Lysis (False): n={len(mech_false)}\n")
        f.write(f"  Median: {mech_false.median():.2f}\n")
        f.write(f"  Mean:   {mech_false.mean():.2f}\n")
        f.write(f"  Std:    {mech_false.std():.2f}\n\n")
        
        f.write(f"Mann-Whitney U Test:\n")
        f.write(f"  Statistic: {stat}\n")
        f.write(f"  P-value:   {p_val:.4e}\n")

    # Plot
    plt.figure(figsize=(8, 6))
    
    # Prepare data for plotting
    plot_data = [mech_true, mech_false]
    labels = ['Mechanical Lysis (True)', 'No Mechanical Lysis (False)']
    colors = ['blue', 'red'] # Blue for True, Red for False
    
    # Create boxplot
    bplot = plt.boxplot(plot_data, labels=labels, patch_artist=True, medianprops=dict(color="black"))
    
    # Color the boxes
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        
    plt.title('Median Reads Remaining (Nonchim) per Dataset by Extraction Type')
    plt.ylabel('Median Number of Reads')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', edgecolor='black', label='Mechanical Lysis: True'),
        Patch(facecolor='red', edgecolor='black', label='Mechanical Lysis: False')
    ]
    plt.legend(handles=legend_elements, loc='best')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "reads_extraction_kit_type.png"))
    plt.close()


def plot_num_asvs_per_dataset(included_datasets_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the Total_Num_ASVs per dataset.
    Exports the results as both a TSV file and a PNG visualization.
    """
    # Filter for datasets with valid Total_Num_ASVs > 0
    df = included_datasets_df.copy()
    df['Total_Num_ASVs'] = pd.to_numeric(df['Total_Num_ASVs'], errors='coerce')
    df = df[df['Total_Num_ASVs'] > 0].dropna(subset=['Total_Num_ASVs'])
    
    # Sort by Total_Num_ASVs descending
    df = df.sort_values('Total_Num_ASVs', ascending=False)
    
    # Prepare data for plotting and export
    plot_data = df[['Dataset_ID', 'Total_Num_ASVs']].set_index('Dataset_ID')
    
    # Optionally export
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        plot_data.to_csv(export_path, sep='\t')
        
    # Create plot
    plt.figure(figsize=(max(10, len(plot_data) * 0.4), 6))
    bars = plt.bar(plot_data.index, plot_data['Total_Num_ASVs'], color='skyblue', edgecolor='black')
    
    plt.title('Total Number of ASVs per Dataset', fontsize=12)
    plt.xlabel('Dataset ID')
    plt.ylabel('Total Number of ASVs')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{int(height)}',
                 ha='center', va='bottom', fontsize=9, rotation=0)
                 
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "num_asvs_per_dataset.png"))
    plt.close()


def plot_num_sgbs_per_dataset(included_datasets_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the Total_Num_SGBs per dataset.
    Exports the results as both a TSV file and a PNG visualization.
    """
    # Filter for datasets with valid Total_Num_SGBs > 0
    df = included_datasets_df.copy()
    df['Total_Num_SGBs'] = pd.to_numeric(df['Total_Num_SGBs'], errors='coerce')
    df = df[df['Total_Num_SGBs'] > 0].dropna(subset=['Total_Num_SGBs'])
    
    # Sort by Total_Num_SGBs descending
    df = df.sort_values('Total_Num_SGBs', ascending=False)
    
    # Prepare data for plotting and export
    plot_data = df[['Dataset_ID', 'Total_Num_SGBs']].set_index('Dataset_ID')
    
    # Optionally export
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        plot_data.to_csv(export_path, sep='\t')
        
    # Create plot
    plt.figure(figsize=(max(10, len(plot_data) * 0.4), 6))
    bars = plt.bar(plot_data.index, plot_data['Total_Num_SGBs'], color='lightgreen', edgecolor='black')
    
    plt.title('Total Number of SGBs per Dataset', fontsize=12)
    plt.xlabel('Dataset ID')
    plt.ylabel('Total Number of SGBs')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{int(height)}',
                 ha='center', va='bottom', fontsize=9, rotation=0)
                 
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "num_sgbs_per_dataset.png"))
    plt.close()


def plot_datasets_by_amplicon_region_coloured_by_technology(included_datasets_df, output_dir, write_tsv_path=None):
    """
    Create a stacked bar plot showing the number of datasets per amplicon region,
    coloured by sequencing technology.
    Replaces NA in Amplicon_Region with "Shotgun" and sorts regions alphabetically.
    """
    # Create a copy and handle NA values
    df = included_datasets_df.copy()
    df['Amplicon_Region'] = df['Amplicon_Region'].fillna('Shotgun')
    
    # Group by Amplicon_Region and Sequencing_Technology
    grouped = df.groupby(['Amplicon_Region', 'Sequencing_Technology']).size().unstack(fill_value=0)
    
    # Sort regions alphabetically
    grouped = grouped.sort_index()
    
    # Optionally export the data used for plotting
    if write_tsv_path is not None:
        export_df = grouped.copy()
        # Rename columns to [Sequencing_Technology]_Count with underscores
        export_df.columns = [f"{col.replace(' ', '_')}_Count" for col in export_df.columns]
        
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        export_df.to_csv(export_path, sep='\t')

    # Create the plot
    plt.figure(figsize=(12, 7))
    
    # Use existing style for plotting
    ax = grouped.plot(kind='bar', stacked=True, figsize=(12, 7), ax=plt.gca())
    
    # Customize the plot
    plt.title('Number of Datasets by Amplicon Region and Sequencing Technology', fontsize=12)
    plt.xlabel('Amplicon Region')
    plt.ylabel('Number of Datasets')
    plt.xticks(rotation=0)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add legend
    plt.legend(title='Sequencing Technology', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, "datasets_by_amplicon_region_coloured_by_technology.png"))
    plt.close()


def plot_age_bucket_proportions(all_samples_df, output_dir, write_tsv_path=None):
    """
    Create a bar plot showing the number of celiac samples (Diagnosed_Celiac=True or Will_Develop_Celiac=True)
    per age bucket. Only includes samples with exact ages (ignores age ranges).
    """
    # Filter for celiac samples (diagnosed or will develop)
    celiac_mask = (all_samples_df['Diagnosed_Celiac'] == True) | (all_samples_df['Will_Develop_Celiac'] == True)
    celiac_samples = all_samples_df[celiac_mask].copy()
    
    # Convert Age to numeric, coercing errors to NaN
    celiac_samples['Age'] = pd.to_numeric(celiac_samples['Age'], errors='coerce')
    
    # Drop samples with NaN ages (these would be age ranges or missing data)
    celiac_samples = celiac_samples.dropna(subset=['Age'])
    
    if celiac_samples.empty:
        print("No celiac samples with exact ages found for age bucket analysis.")
        return
    
    # Initialize bucket counts
    bucket_counts = {bucket_name: 0 for bucket_name in AGE_BUCKETS.keys()}
    
    # Assign each sample to an age bucket
    for age in celiac_samples['Age']:
        for bucket_name, (min_age, max_age) in AGE_BUCKETS.items():
            if min_age <= age <= max_age:
                bucket_counts[bucket_name] += 1
                break
    
    # Create a DataFrame for plotting
    bucket_df = pd.DataFrame({
        'Age_Bucket': list(AGE_BUCKETS.keys()),
        'Num_Samples': [bucket_counts[bucket] for bucket in AGE_BUCKETS.keys()]
    })
    
    # Optionally export the data used for plotting
    if write_tsv_path is not None:
        export_path = write_tsv_path if os.path.isabs(write_tsv_path) else os.path.join(output_dir, write_tsv_path)
        bucket_df.to_csv(export_path, sep='\t', index=False)
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    bars = plt.bar(bucket_df['Age_Bucket'], bucket_df['Num_Samples'], color='skyblue', edgecolor='black')
    
    # Customize the plot
    plt.title('Number of Celiac Samples by Age Bucket', fontsize=12)
    plt.xlabel('Age Bucket')
    plt.ylabel('Number of Samples')
    plt.xticks(rotation=45, ha='right')
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
    plt.savefig(os.path.join(output_dir, "age_bucket_proportions.png"))
    plt.close()
