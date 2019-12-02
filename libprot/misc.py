def safeify_id(record_id: str):
    """Converts all non-alphanumeric characters in `record_id` into hyphens, and returns the string."""
    return ''.join(c if c.isalnum() else '-' for c in record_id)
