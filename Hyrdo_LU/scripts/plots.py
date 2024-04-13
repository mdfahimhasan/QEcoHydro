import numpy as np
import matplotlib.pyplot as plt


def make_Budyko_plot(E_P, PET_P, catchment_no, savepath, E_P_sim=None, PET_P_sim=None):
    """
    Make Budyko plot.

    :param E_P: Evaporative fraction for all catchments.
    :param PET_P: Aridity index for all catchments.
    :param catchment_no: Catchment of interest.
    :param E_P_sim: Evaporative fraction from model simulation.
    :param PET_P_sim: Aridity index from model simulation.
    :param savepath: Plot savepath.

    :return: None.
    """
    # selecting catchment idx
    catchment_idx = catchment_no - 1

    # defining eqn of Budyko curve
    AI_plot_lin = np.arange(0.01, 5.01, 0.01)
    Budyko_curve = lambda AI: (AI * (1 - np.exp(-AI)) * np.tanh(1/AI)) ** 0.5

    # plotting
    plt.plot(PET_P, E_P, 'o', label='catchments')

    plt.plot(PET_P[catchment_idx], E_P[catchment_idx], 'o', label=f'catchment {10} observed')

    if E_P_sim is not None:
        plt.plot(PET_P_sim, E_P_sim, 'ro', label=f'catchment {10} simulated')

    plt.plot(AI_plot_lin, Budyko_curve(AI_plot_lin), color='black', label='Budyko line')

    plt.xlabel('PET/P')
    plt.ylabel('E/P')
    plt.legend()

    plt.savefig(savepath, dpi=100)


def plot_P_PET_Temp(P_mean_monthly_sum, PET_mean_monthly_sum, Temp_mean_monthly, catchment_no, save_path):
    months = list(range(1, 13))
    catchment_idx = catchment_no - 1

    fig, ax = plt.subplots()
    plt.plot(months, P_mean_monthly_sum[:,  catchment_idx], '--g', linewidth=1, label='P')
    plt.plot(months, PET_mean_monthly_sum[:,  catchment_idx], '-g', linewidth=1, label='PET')
    plt.ylabel('mm', color='g')

    ax1 = plt.gca()  # Set the color of the y-axis label and tick labels to red
    ax1.tick_params(axis='y', colors='g')

    ax2 = plt.twinx()  # Create secondary y-axis and plot
    ax2.plot(months, Temp_mean_monthly[:,  catchment_idx], '-r', linewidth=1, label='Temp')
    plt.ylabel('deg. C', color='r')
    ax2.tick_params(axis='y', colors='r')  # Set the color of the secondary y-axis label and tick labels to blue

    handles1, labels1 = ax1.get_legend_handles_labels()  # Save the handles and labels of the first axis
    handles2, labels2 = ax2.get_legend_handles_labels() # Save the handles and labels of the second axis
    handles = handles1 + handles2   # Combine handles and labels from both axes
    labels = labels1 + labels2
    plt.legend(handles, labels, loc='best')
    plt.xlabel('Month')

    fig.savefig(save_path, dpi=100)


def plot_Q_monthly(QT_mean_monthly_observed, QT_mean_monthly_simulated,
                   catchment_no, save_path):
    months = list(range(1, 13))
    catchment_idx = catchment_no - 1

    fig, ax = plt.subplots()
    plt.plot(months, QT_mean_monthly_observed[:, catchment_idx], '-b', linewidth=1, label='Q observed')
    plt.plot(months, QT_mean_monthly_simulated, '-g', linewidth=0.5, label='Q simulated')
    plt.ylabel('mm')
    plt.xlabel('Month')
    plt.legend()

    fig.savefig(save_path, dpi=100)


def plot_Q_daily(QT_daily_observed, QT_daily_simulated, from_day, to_day, year,
                 catchment_no, save_path):

    days = list(range(from_day, to_day))
    catchment_idx = catchment_no - 1

    fig, ax = plt.subplots()
    plt.plot(days, QT_daily_observed[:, catchment_idx][from_day: to_day], '-b', linewidth=0.5, label='Q observed')
    plt.plot(days, QT_daily_simulated[from_day: to_day], '-g', linewidth=0.5, label='Q simulated')

    plt.ylabel('mm')
    plt.xlabel('days')

    tick_labels = list(np.arange(1, to_day + 1 - from_day, 50))
    days = list(range(from_day, to_day, 50))
    plt.xticks(days, tick_labels)
    plt.title(f'for year: {year}')
    plt.legend()

    fig.savefig(save_path, dpi=100)


def plot_simulated_flux_monthly(flux_monthly, color, linestyle, label, ylabel, save_path):
    months = list(range(1, 13))

    fig, ax = plt.subplots()
    plt.plot(months, flux_monthly, color=color, linestyle=linestyle, linewidth=0.5, label=label)
    plt.ylabel(ylabel)
    plt.xlabel('month')
    plt.legend()

    fig.savefig(save_path, dpi=100)


def plot_compare_fluxes(flux_sim_1, flux_sim_2, label1, label2, ylabel, title, save_path):
    months = list(range(1, 13))

    fig, ax = plt.subplots()
    plt.plot(months, flux_sim_1, color='g', linewidth=0.5, label=label1)
    plt.plot(months, flux_sim_2, color='r', linewidth=0.5, label=label2)
    plt.ylabel(ylabel)
    plt.xlabel('month')
    plt.title(title)
    plt.legend()

    fig.savefig(save_path, dpi=100)


def plot_FDC(sim1_QT_daily, sim2_QT_daily, save_path):
    # Sorting flow data in ascending order
    QT_1_sorted = np.sort(sim1_QT_daily)[::-1]  # descending order
    QT_2_sorted = np.sort(sim2_QT_daily)[::-1]  # descending order

    # permanence (similar to exceedance probability)
    permanence = np.arange(1, len(QT_1_sorted) + 1) / len(QT_1_sorted)  # Calculate permanence

    # Semilog plot for flow data
    fig, ax = plt.subplots()
    ax.semilogy(permanence, QT_1_sorted, '-g', linewidth=1, label='SIM 1: 100% Forest')  # Adjust 'color' if needed
    ax.semilogy(permanence, QT_2_sorted, '-r', linewidth=1, label='SIM 2: 70% Forest + 30% grass')  # Adjust 'color' if needed
    # Setting legend and labels
    plt.legend(frameon=False)
    plt.ylabel('Q (mm/d)')
    plt.xlabel('% of Time Exceeded')

    fig.savefig(save_path, dpi=100)


def plot_FFC(dates, sim1_QT_daily, sim2_QT_daily, save_path):
    # Find unique years in the data
    unq_years = np.unique(dates[:, 0])

    # Initialize empty list to store annual maximum flows for both simulations
    max_flow_SIM1 = []
    max_flow_SIM2 = []

    # Iterate over each unique year
    for yr in unq_years:
        # Finding flows of that year
        flows_SIM1 = sim1_QT_daily[dates[:, 0] == yr]
        flows_SIM2 = sim2_QT_daily[dates[:, 0] == yr]

        # Finding maximum flow of that year
        max_SIM1 = max(flows_SIM1)
        max_SIM2 = max(flows_SIM2)

        # Storing to the empty list
        max_flow_SIM1.append(max_SIM1)
        max_flow_SIM2.append(max_SIM2)

    # Sort the annual maximum flows in descending order for both simulations
    max_flow_SIM1 = sorted(max_flow_SIM1, reverse=True)
    max_flow_SIM2 = sorted(max_flow_SIM2, reverse=True)

    # Calculate exceedance probability and return periods
    n_years = len(unq_years)
    exceedance_prob = np.arange(1, n_years + 1) / n_years
    return_periods = 1 / exceedance_prob

    # Plot
    fig, ax = plt.subplots()
    plt.plot(return_periods, max_flow_SIM1, '-g', linewidth=1, label='SIM 1: 100% Forest')
    plt.plot(return_periods, max_flow_SIM2, '-r', linewidth=1, label='SIM 2: 70% Forest + 30% grass')
    plt.legend()
    plt.xlabel('Return Period (years)')
    plt.ylabel('Annual Maximum Streamflow')
    plt.xscale('log')
    plt.yscale('log')

    fig.savefig(save_path, dpi=100)
