"""Constants and enumerations for MDPerTool response time analysis."""

from enum import Enum
from dataclasses import dataclass
from typing import Final


class ModelType(str, Enum):
    """Cumulative fit model types."""
    LOGISTIC = "logistic"
    GOMPERTZ = "gompertz"
    NONE = "none"


class NAReasonCode(str, Enum):
    """Reason codes for N/A (unavailable) fit status."""
    NONE = "none"
    ALL_NONRESPONSIVE = "all_nonresponsive"
    INSUFFICIENT_FIT_POINTS = "insufficient_fit_points"
    MISSING_SIDECAR = "missing_sidecar"


class FitStatus(str, Enum):
    """Fit status indicators."""
    OK = "ok"
    UNAVAILABLE = "unavailable"


@dataclass
class FitConstants:
    """Numerical constants for response time fitting."""
    
    RESPONSE_THRESHOLD: Final[float] = 0.01
    """Energy difference threshold for detecting response."""
    
    CLIP_EXPONENT_MIN: Final[int] = -700
    """Minimum exponent value for numerical stability."""
    
    CLIP_EXPONENT_MAX: Final[int] = 700
    """Maximum exponent value for numerical stability."""
    
    NORMALIZED_EPS: Final[float] = 1e-6
    """Epsilon for normalized values (mask building in curve fitting)."""
    
    SLOPE_TOLERANCE: Final[float] = 1e-10
    """Tolerance for zero slope detection."""
    
    AIC_PARAMETER_COUNT: Final[int] = 2
    """Number of parameters in logistic/Gompertz models."""
