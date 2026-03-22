"""
Residue Response Time Analysis Module

Analyzes per-residue response times from energy decomposition trajectories.
Supports visualization, domain-based clustering, and comparative analysis.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path


class ResidueResponseAnalyzer:
    """Analyze and visualize residue response time dynamics."""

    def __init__(self, response_times_file, metrics_file=None, fit_curve_file=None, 
                 frame_time_delta=1.0):
        """
        Initialize analyzer with response time data.

        Parameters
        ----------
        response_times_file : str
            Path to responseTimes.csv containing per-residue frame indices
        metrics_file : str, optional
            Path to responseTimes_metrics.csv with sigmoid fit parameters
        fit_curve_file : str, optional
            Path to responseTimes_fit_curve.csv with cumulative curves
        frame_time_delta : float
            Time interval per frame (ps), default 1.0
        """
        self.response_times_file = response_times_file
        self.metrics_file = metrics_file
        self.fit_curve_file = fit_curve_file
        self.frame_time_delta = frame_time_delta

        self._load_data()

    def _load_data(self):
        """Load response time data from CSV files."""
        # Load per-residue response times (frame indices)
        response_times_raw = np.loadtxt(self.response_times_file, delimiter=',')
        
        # Handle single residue case
        if response_times_raw.ndim == 0:
            self.response_times_frames = np.array([response_times_raw])
        else:
            self.response_times_frames = response_times_raw.astype(int)

        # Convert frame indices to time (ps)
        self.response_times_ps = self.response_times_frames * self.frame_time_delta
        self.num_residues = len(self.response_times_frames)

        # Load metrics if available
        self.metrics = None
        if self.metrics_file and os.path.isfile(self.metrics_file):
            try:
                self.metrics = pd.read_csv(self.metrics_file).iloc[0].to_dict()
            except Exception:
                pass

        # Load fit curve if available
        self.fit_curve = None
        if self.fit_curve_file and os.path.isfile(self.fit_curve_file):
            try:
                self.fit_curve = pd.read_csv(self.fit_curve_file)
            except Exception:
                pass

    def get_per_residue_summary(self, residue_names=None):
        """
        Get per-residue response time table.

        Parameters
        ----------
        residue_names : list, optional
            Custom residue names/labels. If None, uses residue index.

        Returns
        -------
        pd.DataFrame
            Table with columns: residue_id, residue_name, response_time_frame, response_time_ps
        """
        if residue_names is None:
            residue_names = [f"RES_{i}" for i in range(self.num_residues)]

        data = {
            'residue_id': np.arange(self.num_residues),
            'residue_name': residue_names[:self.num_residues],
            'response_time_frame': self.response_times_frames,
            'response_time_ps': self.response_times_ps,
        }

        return pd.DataFrame(data)

    def get_dynamic_modules(self, time_interval_ps=0.1):
        """
        Classify residues into dynamical modules by response time intervals.

        Parameters
        ----------
        time_interval_ps : float
            Time interval for binning (default 0.1 ps for molecular dynamics)

        Returns
        -------
        dict
            Mapping of module_id -> list of residue indices
        """
        max_time = np.max(self.response_times_ps)
        bins = np.arange(0, max_time + time_interval_ps, time_interval_ps)
        
        modules = {}
        digitized = np.digitize(self.response_times_ps, bins)
        
        for module_id in np.unique(digitized):
            residue_indices = np.where(digitized == module_id)[0].tolist()
            if residue_indices:
                modules[int(module_id)] = residue_indices

        return modules

    def get_domain_summary(self, domain_residue_mapping):
        """
        Summarize response times by structural domain.

        Parameters
        ----------
        domain_residue_mapping : dict
            Mapping of domain_name -> list of residue indices

        Returns
        -------
        pd.DataFrame
            Domain summary with statistics
        """
        results = []
        
        for domain_name, residue_indices in domain_residue_mapping.items():
            valid_indices = [i for i in residue_indices if i < self.num_residues]
            
            if not valid_indices:
                continue

            times_ps = self.response_times_ps[valid_indices]
            times_frames = self.response_times_frames[valid_indices]

            results.append({
                'domain': domain_name,
                'num_residues': len(valid_indices),
                'min_response_time_ps': float(np.min(times_ps)),
                'max_response_time_ps': float(np.max(times_ps)),
                'mean_response_time_ps': float(np.mean(times_ps)),
                'median_response_time_ps': float(np.median(times_ps)),
                'std_response_time_ps': float(np.std(times_ps)),
                'min_response_frame': int(np.min(times_frames)),
                'max_response_frame': int(np.max(times_frames)),
            })

        return pd.DataFrame(results)

    def get_sigmoid_parameters(self):
        """
        Get sigmoid curve fitting parameters (if available).

        Returns
        -------
        dict
            Parameters: t_half_empirical, t_half_fit, k_d (slope), rmse, model type
        """
        if self.metrics is None:
            return None

        return {
            'half_response_time_empirical': self.metrics.get('t_half_empirical_frame'),
            'half_response_time_fit': self.metrics.get('t_half_fit_frame'),
            'dissipation_rate_constant': self.metrics.get('k_d'),
            'fit_rmse': self.metrics.get('fit_rmse'),
            'model_type': self.metrics.get('selected_model'),
            'total_residues': self.metrics.get('total_residues'),
            'responded_residues': self.metrics.get('responded_residues'),
            'responded_fraction': self.metrics.get('responded_fraction'),
        }

    def get_cumulative_curve_data(self):
        """
        Get cumulative response curve and fitted models.

        Returns
        -------
        dict or None
            Keys: 'frames', 'cumulative_observed', 'cumulative_logistic', 
                  'cumulative_gompertz', 'cumulative_fitted', 'model_type'
        """
        if self.fit_curve is None:
            return None

        return {
            'frames': self.fit_curve['frame'].values,
            'cumulative_observed': self.fit_curve['cumulative_responded_observed'].values,
            'cumulative_logistic': self.fit_curve['cumulative_responded_logistic'].values,
            'cumulative_gompertz': self.fit_curve['cumulative_responded_gompertz'].values,
            'cumulative_fitted': self.fit_curve['cumulative_responded_fitted'].values,
            'model_type': self.fit_curve['selected_model'].iloc[0] if len(self.fit_curve) > 0 else None,
        }

    def get_residue_groups_by_threshold(self, fast_threshold_ps=None, slow_threshold_ps=None):
        """
        Classify residues by response speed.

        Parameters
        ----------
        fast_threshold_ps : float, optional
            Upper threshold for "fast" response (default: median response time)
        slow_threshold_ps : float, optional
            Lower threshold for "slow" response (default: median response time)

        Returns
        -------
        dict
            Keys: 'fast', 'medium', 'slow' with residue indices
        """
        if fast_threshold_ps is None:
            fast_threshold_ps = np.median(self.response_times_ps) * 0.5

        if slow_threshold_ps is None:
            slow_threshold_ps = np.median(self.response_times_ps) * 1.5

        fast_residues = np.where(self.response_times_ps <= fast_threshold_ps)[0].tolist()
        slow_residues = np.where(self.response_times_ps >= slow_threshold_ps)[0].tolist()
        medium_residues = [i for i in range(self.num_residues) 
                          if i not in fast_residues and i not in slow_residues]

        return {
            'fast': fast_residues,
            'medium': medium_residues,
            'slow': slow_residues,
            'fast_threshold_ps': fast_threshold_ps,
            'slow_threshold_ps': slow_threshold_ps,
        }

    def export_analysis_summary(self, output_dir, residue_names=None):
        """
        Export comprehensive analysis summary to CSV files.

        Parameters
        ----------
        output_dir : str
            Output directory for CSV files
        residue_names : list, optional
            Custom residue names
        """
        os.makedirs(output_dir, exist_ok=True)

        # Per-residue table
        per_residue = self.get_per_residue_summary(residue_names)
        per_residue.to_csv(os.path.join(output_dir, 'residue_response_summary.csv'), index=False)

        # Response groups
        groups = self.get_residue_groups_by_threshold()
        groups_df = pd.DataFrame({
            'category': ['fast', 'medium', 'slow'],
            'count': [len(groups['fast']), len(groups['medium']), len(groups['slow'])],
            'threshold_ps': [f"<= {groups['fast_threshold_ps']:.2f}", 
                            f"> {groups['fast_threshold_ps']:.2f} and < {groups['slow_threshold_ps']:.2f}",
                            f">= {groups['slow_threshold_ps']:.2f}"],
        })
        groups_df.to_csv(os.path.join(output_dir, 'residue_response_groups.csv'), index=False)

        # Sigmoid parameters
        if self.metrics:
            metrics_df = pd.DataFrame([self.metrics])
            metrics_df.to_csv(os.path.join(output_dir, 'sigmoid_parameters.csv'), index=False)

        return True


def load_multiple_analyses(analysis_output_dirs, frame_time_delta=1.0):
    """
    Load and compare response times from multiple analyses.

    Parameters
    ----------
    analysis_output_dirs : list
        List of analysis output directories containing responseTimes.csv files
    frame_time_delta : float
        Time interval per frame (ps)

    Returns
    -------
    dict
        Mapping of analysis_name -> ResidueResponseAnalyzer
    """
    analyzers = {}
    
    for analysis_dir in analysis_output_dirs:
        response_file = os.path.join(analysis_dir, 'responseTimes.csv')
        metrics_file = os.path.join(analysis_dir, 'responseTimes_metrics.csv')
        fit_file = os.path.join(analysis_dir, 'responseTimes_fit_curve.csv')

        if os.path.isfile(response_file):
            analysis_name = os.path.basename(analysis_dir)
            try:
                analyzer = ResidueResponseAnalyzer(
                    response_file, 
                    metrics_file if os.path.isfile(metrics_file) else None,
                    fit_file if os.path.isfile(fit_file) else None,
                    frame_time_delta
                )
                analyzers[analysis_name] = analyzer
            except Exception as e:
                print(f"Error loading analysis {analysis_name}: {e}")

    return analyzers


if __name__ == '__main__':
    # Example usage
    analyzer = ResidueResponseAnalyzer(
        'responseTimes.csv',
        'responseTimes_metrics.csv',
        'responseTimes_fit_curve.csv',
        frame_time_delta=1.0
    )

    print("Per-residue summary:")
    print(analyzer.get_per_residue_summary())

    print("\nSigmoid parameters:")
    print(analyzer.get_sigmoid_parameters())

    print("\nDynamic modules:")
    print(analyzer.get_dynamic_modules())
