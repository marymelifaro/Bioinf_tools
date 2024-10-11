from numbers import Real

from config import PHRED33_SHIFT


def is_bounded(x, bounds):
    return bounds[0] <= x <= bounds[1]


def filter_length(seq: str, length_bounds: tuple[int, int]) -> bool:
    return is_bounded(len(seq), length_bounds)


def filter_gc(seq: str, gc_bounds: tuple[Real, Real]) -> bool:
    gc_value = (seq.count('G') + seq.count('C')) / len(seq) * 100
    return is_bounded(gc_value, gc_bounds)


def filter_quality(quality: str, quality_threshold: Real) -> bool:
    quality_value = sum(ord(x) - PHRED33_SHIFT for x in quality) / len(quality)
    return quality_value >= quality_threshold
