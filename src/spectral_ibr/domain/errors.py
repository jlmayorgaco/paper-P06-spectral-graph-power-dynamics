"""Domain-specific exceptions."""


class SpectralIbrError(Exception):
    """Base exception for domain failures."""


class InvalidLinearizedModelError(SpectralIbrError):
    """Raised when a partitioned state matrix is dimensionally invalid."""


class IbrConversionError(SpectralIbrError):
    """Raised when an SG-to-IBR conversion cannot be completed."""


class MissingBaseCaseError(IbrConversionError):
    """Raised when the configured base workbook is not present."""


class InvalidConversionProfileError(IbrConversionError):
    """Raised when the requested IBR conversion profile is invalid."""
