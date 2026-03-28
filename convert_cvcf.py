import os, json
base = r'c:\Users\admin\codes\CEWorkbench\data\CEWorkbench\hamiltonians'
systems = ['Nb-Ti','Nb-V','Nb-Zr','Ti-V','Ti-Zr','V-Zr']
for sys in systems:
    src = os.path.join(base, f'{sys}_BCC_A2_T', 'hamiltonian.json')
    if not os.path.exists(src):
        print('source missing', src)
        continue
    data = json.load(open(src, 'r', encoding='utf-8'))
    cft = {t['name']: t for t in data.get('cecTerms', [])}
    v0 = float(cft.get('CF_0', {}).get('a', 0.0))
    v1 = float(cft.get('CF_1', {}).get('a', 0.0))
    v2 = float(cft.get('CF_2', {}).get('a', 0.0))
    v3 = float(cft.get('CF_3', {}).get('a', 0.0))
    cvcf = {
        'elements': data.get('elements', ''),
        'structurePhase': data.get('structurePhase', 'BCC_A2'),
        'model': 'CVCF',
        'cecTerms': [
            {'name': 'e4AB',  'description': '4-site tetrahedron CF', 'numSites': 4, 'multiplicity': 6.0, 'a': -16.0*v0, 'b': 0.0},
            {'name': 'e3AB',  'description': '3-site triangle CF', 'numSites': 3, 'multiplicity': 12.0, 'a': v1, 'b': 0.0},
            {'name': 'e22AB', 'description': '2-site pair CF (4 bonds)', 'numSites': 2, 'multiplicity': 4.0, 'a': -16.0*v2, 'b': 0.0},
            {'name': 'e21AB', 'description': '2-site pair CF (3 bonds)', 'numSites': 2, 'multiplicity': 3.0, 'a': -12.0*v3, 'b': 0.0}
        ],
        'cecUnits': data.get('cecUnits', 'J/mol'),
        'reference': data.get('reference', ''),
        'notes': 'Auto-converted to CVCF basis from T orthogonal data',
        'ncf': 4
    }
    dst_dir = os.path.join(base, f'{sys}_BCC_A2_T_CVCF')
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, 'hamiltonian.json')
    with open(dst, 'w', encoding='utf-8') as f:
        json.dump(cvcf, f, indent=2)
    print('wrote', dst)
