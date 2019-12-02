def safeify_id(record_id: str):
    return ''.join(c if c.isalnum() else '-' for c in record_id)
