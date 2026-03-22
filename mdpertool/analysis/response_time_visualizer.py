"""
Response Time Visualization Helper

Provides matplotlib-based visualization for residue response time analysis,
including sigmoid curves, cumulative histograms, and per-residue distributions.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd


class ResponseTimeVisualizer:
    """Visualize response time dynamics with various plot types."""

    @staticmethod
    def plot_cumulative_response_curve(ax, cumulative_data, title="Energy Dissipation Curve"):
        """
        Plot cumulative response curve with fitted models.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot on
        cumulative_data : dict
            Output from ResidueResponseAnalyzer.get_cumulative_curve_data()
            Keys: 'frames', 'cumulative_observed', 'cumulative_logistic', 
                  'cumulative_gompertz', 'cumulative_fitted', 'model_type'
        title : str
            Plot title
        """
        if cumulative_data is None:
            ax.text(0.5, 0.5, 'No cumulative curve data available', 
                   ha='center', va='center', transform=ax.transAxes)
            return

        frames = cumulative_data['frames']
        cumulative_obs = cumulative_data['cumulative_observed']
        cumulative_logistic = cumulative_data['cumulative_logistic']
        cumulative_gompertz = cumulative_data['cumulative_gompertz']
        cumulative_fitted = cumulative_data['cumulative_fitted']
        model_type = cumulative_data['model_type']

        # Plot observed data
        ax.plot(frames, cumulative_obs, 'o-', color='black', linewidth=2, 
                markersize=4, label='Observed', zorder=5)

        # Plot fitted models
        valid_logistic = ~np.isnan(cumulative_logistic)
        valid_gompertz = ~np.isnan(cumulative_gompertz)

        if np.any(valid_logistic):
            ax.plot(frames[valid_logistic], cumulative_logistic[valid_logistic], 
                   '--', color='blue', linewidth=1.5, label='Logistic Fit', alpha=0.7)

        if np.any(valid_gompertz):
            ax.plot(frames[valid_gompertz], cumulative_gompertz[valid_gompertz], 
                   '--', color='green', linewidth=1.5, label='Gompertz Fit', alpha=0.7)

        # Highlight selected model
        if model_type and model_type != 'none':
            if model_type == 'logistic':
                valid_fitted = ~np.isnan(cumulative_fitted)
                if np.any(valid_fitted):
                    ax.plot(frames[valid_fitted], cumulative_fitted[valid_fitted], 
                           '-', color='red', linewidth=2.5, label='Selected Model', zorder=4)
            elif model_type == 'gompertz':
                valid_fitted = ~np.isnan(cumulative_fitted)
                if np.any(valid_fitted):
                    ax.plot(frames[valid_fitted], cumulative_fitted[valid_fitted], 
                           '-', color='red', linewidth=2.5, label='Selected Model', zorder=4)

        ax.set_xlabel('Time Frame', fontsize=11)
        ax.set_ylabel('Cumulative # Residues Responded', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(loc='lower right', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    @staticmethod
    def plot_response_time_histogram(ax, response_times_ps, title="Per-Residue Response Time Distribution"):
        """
        Plot histogram of response times across residues.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot on
        response_times_ps : array-like
            Per-residue response times in picoseconds
        title : str
            Plot title
        """
        if response_times_ps is None or len(response_times_ps) == 0:
            ax.text(0.5, 0.5, 'No response time data available', 
                   ha='center', va='center', transform=ax.transAxes)
            return

        # Calculate histogram
        bins = max(10, len(response_times_ps) // 5)  # Adaptive bin count
        counts, edges, patches = ax.hist(response_times_ps, bins=bins, 
                                        color='steelblue', alpha=0.7, edgecolor='black')

        # Color code by quartiles
        q1 = np.percentile(response_times_ps, 25)
        q2 = np.percentile(response_times_ps, 50)
        q3 = np.percentile(response_times_ps, 75)

        for i, patch in enumerate(patches):
            center = (edges[i] + edges[i+1]) / 2.0
            if center <= q1:
                patch.set_facecolor('#77dd77')  # Green: fast
            elif center <= q2:
                patch.set_facecolor('#77dd77')
            elif center <= q3:
                patch.set_facecolor('#f6bc66')  # Orange: medium
            else:
                patch.set_facecolor('#ff686b')  # Red: slow

        # Add vertical lines for quartiles
        ax.axvline(q1, color='green', linestyle='--', linewidth=1, alpha=0.7, label=f'Q1: {q1:.1f} ps')
        ax.axvline(q2, color='orange', linestyle='--', linewidth=1, alpha=0.7, label=f'Q2: {q2:.1f} ps')
        ax.axvline(q3, color='red', linestyle='--', linewidth=1, alpha=0.7, label=f'Q3: {q3:.1f} ps')

        ax.set_xlabel('Response Time (ps)', fontsize=11)
        ax.set_ylabel('# Residues', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(True, alpha=0.3, axis='y')

    @staticmethod
    def plot_response_timeline(ax, response_times_ps, residue_ids=None, 
                              title="Residue Response Timeline"):
        """
        Plot response times as a timeline, colored by speed.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot on
        response_times_ps : array-like
            Per-residue response times
        residue_ids : array-like, optional
            Custom residue identifiers
        title : str
            Plot title
        """
        if response_times_ps is None or len(response_times_ps) == 0:
            ax.text(0.5, 0.5, 'No response time data available', 
                   ha='center', va='center', transform=ax.transAxes)
            return

        if residue_ids is None:
            residue_ids = np.arange(len(response_times_ps))

        # Sort by response time for visual clarity
        sorted_indices = np.argsort(response_times_ps)
        sorted_times = response_times_ps[sorted_indices]
        sorted_ids = np.array(residue_ids)[sorted_indices]

        # Color by quartile
        q1 = np.percentile(response_times_ps, 25)
        q3 = np.percentile(response_times_ps, 75)

        colors = []
        for t in sorted_times:
            if t <= q1:
                colors.append('#77dd77')  # Fast: green
            elif t >= q3:
                colors.append('#ff686b')  # Slow: red
            else:
                colors.append('#f6bc66')  # Medium: orange

        # Plot as horizontal bars
        y_pos = np.arange(len(sorted_ids))
        ax.barh(y_pos, sorted_times, color=colors, edgecolor='black', linewidth=0.5)

        # Sample residue labels for clarity (show every Nth)
        step = max(1, len(sorted_ids) // 20)
        tick_positions = y_pos[::step]
        tick_labels = [str(sorted_ids[i]) for i in range(0, len(sorted_ids), step)]

        ax.set_yticks(tick_positions)
        ax.set_yticklabels(tick_labels, fontsize=8)
        ax.set_xlabel('Response Time (ps)', fontsize=11)
        ax.set_ylabel('Residue ID (sorted)', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')

    @staticmethod
    def plot_response_modules(ax, response_times_frames, module_interval=1, 
                             title="Dynamic Modules by Response Time"):
        """
        Plot residues grouped into dynamical modules.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot on
        response_times_frames : array-like
            Per-residue response times in frame units
        module_interval : float
            Frame interval for each module
        title : str
            Plot title
        """
        if response_times_frames is None or len(response_times_frames) == 0:
            ax.text(0.5, 0.5, 'No response data available', 
                   ha='center', va='center', transform=ax.transAxes)
            return

        # Bin residues into modules
        max_frame = np.max(response_times_frames)
        num_modules = int(np.ceil(max_frame / module_interval)) + 1

        module_counts = np.zeros(num_modules)
        for frame_idx in response_times_frames:
            module_idx = int(frame_idx / module_interval)
            if module_idx < num_modules:
                module_counts[module_idx] += 1

        # Plot module distribution
        module_times = np.arange(num_modules) * module_interval
        colors = plt.cm.viridis(np.linspace(0, 1, num_modules))

        bars = ax.bar(module_times, module_counts, width=module_interval * 0.8, 
                     color=colors, edgecolor='black', linewidth=1)

        ax.set_xlabel('Module Response Time (frames)', fontsize=11)
        ax.set_ylabel('# Residues', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_xlim(left=-module_interval/2)


if __name__ == '__main__':
    # Example usage
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Generate sample data
    np.random.seed(42)
    response_times = np.random.gamma(shape=2, scale=100, size=200) + 50
    frames = (response_times / 1.0).astype(int)

    # Sample cumulative data
    frames_arr = np.arange(0, np.max(frames) + 100, 10)
    cumulative_sample = {
        'frames': frames_arr,
        'cumulative_observed': np.sqrt(frames_arr) * 2,
        'cumulative_logistic': np.sqrt(frames_arr) * 2.1,
        'cumulative_gompertz': np.sqrt(frames_arr) * 2.05,
        'cumulative_fitted': np.sqrt(frames_arr) * 2.05,
        'model_type': 'gompertz'
    }

    viz = ResponseTimeVisualizer()
    viz.plot_cumulative_response_curve(axes[0, 0], cumulative_sample)
    viz.plot_response_time_histogram(axes[0, 1], response_times)
    viz.plot_response_timeline(axes[1, 0], response_times)
    viz.plot_response_modules(axes[1, 1], frames)

    plt.tight_layout()
    plt.show()
