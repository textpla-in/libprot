def safeify_id(record_id: str):
    return ''.join(c if c.isalnum() else '_' for c in record_id)
