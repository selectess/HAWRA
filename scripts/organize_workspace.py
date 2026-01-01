import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

def iter_dirs(root):
    for base, dirs, files in os.walk(root):
        # skip virtualenv and dist zips
        if os.path.basename(base) in {'.venv'}:
            continue
        yield base, dirs, files

def main():
    created = []
    for base, dirs, files in iter_dirs(ROOT):
        if base.endswith('__pycache__'):
            continue
        # consider empty if no files except hidden and no subdirs
        visible_files = [f for f in files if not f.startswith('.')]
        if not visible_files and len(dirs) == 0:
            keep_path = os.path.join(base, '.keep')
            if not os.path.exists(keep_path):
                with open(keep_path, 'w') as f:
                    f.write('')
                created.append(keep_path)
    print('KEEP created:', len(created))

if __name__ == '__main__':
    main()