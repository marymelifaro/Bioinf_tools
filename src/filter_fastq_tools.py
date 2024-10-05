PHRED33_SHIFT = 33
number = int | float


def filter_length(seq: str, length_bounds: tuple[int, int]) -> bool:
    return len(seq) < length_bounds[0] or len(seq) > length_bounds[1]


def filter_gc(seq: str, gc_bounds: tuple[number, number]) -> bool:
    gc_value = (seq.count('G') + seq.count('C')) / len(seq) * 100
    return gc_value < gc_bounds[0] or gc_value > gc_bounds[1]


def filter_quality(quality: str, quality_threshold: number) -> bool:
    quality_value = sum(ord(x) - PHRED33_SHIFT for x in quality) / len(quality)
    return quality_value < quality_threshold
