#print('Hello World!')

def get_data_path():
    import os
    path = os.path.join(os.path.dirname(__file__), 'data')
    return path
    