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

    plt.plot(AI_plot_lin, Budyko_curve(AI_plot_lin), label='Budyko line')

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
    ax2.plot(months, Temp_mean_monthly[:,  catchment_idx], '-b', linewidth=1, label='Temp')
    plt.ylabel('deg. C', color='b')
    ax2.tick_params(axis='y', colors='b')  # Set the color of the secondary y-axis label and tick labels to blue

    handles1, labels1 = ax1.get_legend_handles_labels()  # Save the handles and labels of the first axis
    handles2, labels2 = ax2.get_legend_handles_labels() # Save the handles and labels of the second axis
    handles = handles1 + handles2   # Combine handles and labels from both axes
    labels = labels1 + labels2
    plt.legend(handles, labels, loc='best')
    plt.xlabel('Month')

    fig.savefig(save_path, dpi=100)


def plot_Q_monthly(QT_mean_monthly_observed, QT_mean_monthly_simulated,
                   QS_mean_monthly_simulated, QF_mean_monthly_simulated,
                   catchment_no, save_path):
    months = list(range(1, 13))
    catchment_idx = catchment_no - 1

    fig, ax = plt.subplots()
    plt.plot(months, QT_mean_monthly_observed[:, catchment_idx], '-g', linewidth=1, label='Q observed')
    plt.plot(months, QT_mean_monthly_simulated, linewidth=0.5, label='Q simulated')
    plt.plot(months, QS_mean_monthly_simulated, '--b', linewidth=0.5, label='Qs simulated')
    plt.plot(months, QF_mean_monthly_simulated, '--r', linewidth=0.5, label='Qf simulated')
    plt.ylabel('mm')
    plt.xlabel('Month')
    plt.legend()

    fig.savefig(save_path, dpi=100)


def plot_Q_daily(QT_daily_observed, QT_daily_simulated, from_day, to_day, year,
                 catchment_no, save_path):

    days = list(range(from_day, to_day))
    catchment_idx = catchment_no - 1

    fig, ax = plt.subplots()
    plt.plot(days, QT_daily_observed[:, catchment_idx][from_day: to_day], '-g', linewidth=0.5, label='Q observed')
    plt.plot(days, QT_daily_simulated[from_day: to_day], linewidth=0.5, label='Q simulated')

    plt.ylabel('mm')
    plt.xlabel('days')

    tick_labels = list(np.arange(1, to_day + 1 - from_day, 50))
    days = list(range(from_day, to_day, 50))
    plt.xticks(days, tick_labels)
    plt.title(f'for year: {year}')
    plt.legend()

    fig.savefig(save_path, dpi=100)


def plot_Q_SM_daily(SM_daily_simulated, QT_daily_simulated, from_day, to_day, year,
                    save_path):

    days = list(range(from_day, to_day))

    fig, ax = plt.subplots()
    plt.plot(days, SM_daily_simulated[from_day: to_day], '-g', linewidth=0.5, label='Soil Moisture')
    ax1 = plt.gca()  # Set the color of the y-axis label and tick labels to red
    ax1.tick_params(axis='y', colors='g')
    plt.ylabel('mm', color='g')
    plt.xlabel('days')

    ax2 = plt.twinx()  # Create secondary y-axis and plot
    ax2.plot(days, QT_daily_simulated[from_day: to_day], linewidth=0.5, label='Q simulated')
    plt.ylabel('mm', color='b')
    ax2.tick_params(axis='y', colors='b')  # Set the color of the secondary y-axis label and tick labels to blue

    handles1, labels1 = ax1.get_legend_handles_labels()  # Save the handles and labels of the first axis
    handles2, labels2 = ax2.get_legend_handles_labels() # Save the handles and labels of the second axis
    handles = handles1 + handles2   # Combine handles and labels from both axes
    labels = labels1 + labels2

    tick_labels = list(np.arange(1, to_day + 1 - from_day, 50))
    days = list(range(from_day, to_day, 50))
    plt.xticks(days, tick_labels)
    plt.title(f'for year: {year}')
    plt.legend(handles, labels, loc='best')

    fig.savefig(save_path, dpi=100)