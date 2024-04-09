import numpy as np
import pandas as pd


def make_means_new(input, dates):
    """
    Aggregates flux values as mean of monthly total for all years and annual mean.
    """
    years = np.unique(dates[:, 0])
    years = np.delete(years, 0)  # removing 1st years (probably because its warm-up period)

    # units in annual mean (something/year)
    input_annual = np.array([np.nansum(input[dates[:, 0] == year], axis=0) for year in years])
    input_mean_annual = np.nanmean(input_annual, axis=0)

    # units in mean of monthly total (something/month) over all years
    input_mean_monthly = np.array([np.nansum(input[dates[:, 1] == month], axis=0) / len(years) for month in range(1, 13)])

    return input_mean_monthly, input_mean_annual


def make_means(input, dates):
    """
    Aggregates flux values as daily mean across a month and annual mean.
    """

    years = np.unique(dates[:, 0])
    years = np.delete(years, 0)
    input_annual = []

    # units in annual mean (something/year)
    for year in years:
        aux = np.where(dates[:, 0] == year)[0]
        #print(f"Indices for year {year}: {aux}")  # Check the indices being used for annual sum
        # Check if input is 1D or 2D
        if input.ndim == 1:
            input_annual.append(np.nansum(input[aux], axis=0))  # For 1D input
        else:
            input_annual.append(np.nansum(input[aux, :], axis=0))  # For 2D input

    input_annual = np.array(input_annual)
    input_mean_annual = np.nanmean(input_annual, axis=0)

    # units in monthly mean (something/day)
    input_mean_monthly = []
    for i in range(1, 13):  # Adjusted for 0-based indexing
        aux = np.where(dates[:, 1] == i)[0]

        # Check if input is 1D or 2D
        if input.ndim == 1:
            input_mean_monthly.append(np.nanmean(input[aux], axis=0))  # For 1D input
        else:
            input_mean_monthly.append(np.nanmean(input[aux, :], axis=0))  # For 2D input

    input_mean_monthly = np.array(input_mean_monthly)

    return input_mean_monthly, input_mean_annual


def make_monthly_new(input, dates, option):
    """
    Aggregates flux values as mean (something/day) or sum (something/month) of flux values
    for a specific year and a month.

    :param input: The input data comes from model() function in model.py.
                  It is a list of daily values across all years and months
    :param dates: A numpy array consisting of years, months, and dates information.
    :param option: str indicating what aggregation operations to perform. Can be one of 'mean' or 'sum'.
    """
    # converting dates to a pandas DataFrame for easier manipulation
    df = pd.DataFrame(dates, columns=['year', 'month', 'day'])

    # the input data comes from model() function in model.py
    # the input is a list of daily values across all years and months
    # adding the input list as a column in the pandas dataframe
    df['flux'] = pd.Series(input)

    # years and months
    unq_years = df['year'].unique()
    months = list(range(1, 13))

    monthly_results = []  # an empty list to store results of monthly aggregation

    for year in unq_years:
        for month in months:
            temp_df = df[(df.year == year) & (df.month == month)]

            if option == 'mean':
                result = np.nanmean(temp_df['flux'])  # unit in something/day

            elif option == 'sum':
                result = np.nansum(temp_df['flux'])  # unit in something/month

            monthly_results.append(result)

    monthly_results = np.vstack(monthly_results)

    return monthly_results


def aggregate(output_dict, dates):
    # Aggregate Monthly
    output_dict['Ea_mon'] = make_monthly_new(output_dict['Ea'], dates, option='sum')
    output_dict['QF_mon'] = make_monthly_new(output_dict['QF'], dates, option='sum')
    output_dict['R_mon'] = make_monthly_new(output_dict['R'], dates, option='sum')
    output_dict['QS_mon'] = make_monthly_new(output_dict['QS'], dates, option='sum')
    output_dict['QT_mon'] = make_monthly_new(output_dict['QT'], dates, option='sum')
    output_dict['Su_mon'] = make_monthly_new(output_dict['Su'], dates, option='mean')
    output_dict['Ss_mon'] = make_monthly_new(output_dict['Ss'], dates, option='mean')
    output_dict['St_mon'] = make_monthly_new(output_dict['St'], dates, option='mean')
    output_dict['IE_mon'] = make_monthly_new(output_dict['IE'], dates, option='sum')
    output_dict['SE_mon'] = make_monthly_new(output_dict['SE'], dates, option='sum')
    output_dict['P_mon'] = make_monthly_new(output_dict['P'], dates, option='sum')
    output_dict['Ei_mon'] = make_monthly_new(output_dict['Ei'], dates, option='sum')
    output_dict['Et_mon'] = make_monthly_new(output_dict['Et'], dates, option='sum')
    output_dict['PET_mon'] = make_monthly_new(output_dict['Ep'], dates, option='sum')
    output_dict['S_canopy_mon'] = make_monthly_new(output_dict['S_canopy'], dates, option='mean')
    output_dict['pot_inf_mon'] = make_monthly_new(output_dict['pot_inf'], dates, option='sum')

    # Aggregate mean_monthly and mean annual
    output_dict['E_mean_mon'], output_dict['E_mean_annual'] = make_means_new(output_dict['Ea'], dates)
    output_dict['QF_mean_mon'], output_dict['QF_mean_annual'] = make_means_new(output_dict['QF'], dates)
    output_dict['R_mean_mon'], output_dict['R_mean_annual'] = make_means_new(output_dict['R'], dates)
    output_dict['QS_mean_mon'], output_dict['QS_mean_annual'] = make_means_new(output_dict['QS'], dates)
    output_dict['QT_mean_mon'], output_dict['QT_mean_annual'] = make_means_new(output_dict['QT'], dates)
    output_dict['Su_mean_mon'], output_dict['Su_mean_annual'] = make_means(output_dict['Su'], dates)
    output_dict['Ss_mean_mon'], output_dict['Ss_mean_annual'] = make_means(output_dict['Ss'], dates)
    output_dict['St_mean_mon'], output_dict['St_mean_annual'] = make_means(output_dict['St'], dates)
    output_dict['IE_mean_mon'], output_dict['IE_mean_annual'] = make_means_new(output_dict['IE'], dates)
    output_dict['SE_mean_mon'], output_dict['SE_mean_annual'] = make_means_new(output_dict['SE'], dates)
    output_dict['P_mean_mon'], output_dict['P_mean_annual'] = make_means_new(output_dict['P'], dates)
    output_dict['Ei_mean_mon'], output_dict['Ei_mean_annual'] = make_means_new(output_dict['Ei'], dates)
    output_dict['Et_mean_mon'], output_dict['Et_mean_annual'] = make_means_new(output_dict['Et'], dates)
    output_dict['PET_mean_mon'], output_dict['PET_mean_annual'] = make_means_new(output_dict['Ep'], dates)
    output_dict['S_canopy_mean_mon'], output_dict['S_canopy_mean_annual'] = make_means_new(output_dict['S_canopy'], dates)
    output_dict['pot_inf_mean_mon'], output_dict['pot_inf_mean_annual'] = make_means_new(output_dict['pot_inf'], dates)

    return output_dict