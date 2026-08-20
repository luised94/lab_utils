import pandas as pd, numpy as np, tifffile, json, base64, io
from PIL import Image

prof = pd.read_csv('/mnt/user-data/uploads/manual_lane_profiles.csv')
band = pd.read_csv('/mnt/user-data/uploads/band_measurements.csv')
img  = tifffile.imread('/mnt/user-data/uploads/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif')

N = 1000
roi_y = int(prof.roi_y.iloc[0]); roi_h = int(prof.roi_h.iloc[0])
mm_per_px = float((prof.groupby('lane_index').migration_position_millimetres.last().iloc[0]) / (N-1))
bg = float(prof.plate_background_median.iloc[0])

# per-lane profile (raw_value), roi_x/roi_w
lanes = sorted(prof.lane_index.unique())
profiles = {}
roi = {}
for li in lanes:
    s = prof[prof.lane_index==li].sort_values('migration_position_pixels')
    profiles[str(li)] = [int(round(v)) for v in s.raw_value.values]
    roi[str(li)] = {'x': int(s.roi_x.iloc[0]), 'w': int(s.roi_w.iloc[0])}

# identities (measured lanes only, from band_measurements)
ident = {}
for li,sub in band.groupby('lane_index'):
    r=sub.iloc[0]
    ident[str(int(li))] = {'well': int(r.well_number), 'label': r.sample_label, 'role': r.analysis_role}

# consensus band metadata + per-lane per-band values
bands_meta = []
for bi,sub in band.groupby('canonical_band_index'):
    r0 = sub.iloc[0]
    fragile = int((sub.baseline_agreement_status=='fragile').sum())
    bands_meta.append({
        'index': int(bi),
        'mm': float(r0.canonical_position_millimetres),
        'ws': int(r0.window_start_pixels),
        'we': int(r0.window_end_pixels),
        'occupancy': int(r0.band_occupancy),
        'supported': (r0.band_is_well_supported=='yes'),
        'spread_mm': float(r0.cross_lane_spread_millimetres),
        'fragile_lanes': fragile,
    })
bands_meta.sort(key=lambda d: d['index'])

# per (lane, band): area, apex, fragile flag, apex migration px
per = {}
for _,r in band.iterrows():
    per.setdefault(str(int(r.lane_index)), {})[str(int(r.canonical_band_index))] = {
        'area': float(r.reported_area),
        'apex': float(r.apex_height_above_baseline),
        'fragile': (r.baseline_agreement_status=='fragile'),
        'apex_px': int(r.apex_migration_pixels),
        'basis': r.reported_value_basis,
    }

# gel crop image: rows roi_y..roi_y+roi_h, cols spanning all lane ROIs
x0 = min(roi[str(li)]['x'] for li in lanes)
x1 = max(roi[str(li)]['x']+roi[str(li)]['w'] for li in lanes)
crop = img[roi_y:roi_y+roi_h, x0:x1].astype(np.float32)
vmin = np.percentile(img,1.0); vmax = np.percentile(img,99.5)
disp = np.clip((crop - vmin)/(vmax-vmin), 0, 1)
gray_r = (255*(1.0-disp)).astype(np.uint8)   # dark bands on white
pil = Image.fromarray(gray_r, mode='L')
# downscale width for payload size, keep height resolution reasonable
target_w = 480
target_h = int(round(pil.height * target_w/pil.width))
pil = pil.resize((target_w, target_h), Image.BILINEAR)
buf = io.BytesIO(); pil.save(buf, format='PNG', optimize=True)
gel_b64 = base64.b64encode(buf.getvalue()).decode('ascii')

payload = {
    'gel_id': band.gel_id.iloc[0],
    'N': N, 'mm_per_px': mm_per_px, 'roi_y': roi_y, 'bg': bg,
    'lanes': [int(x) for x in lanes],
    'x0': int(x0), 'x1': int(x1),
    'profiles': profiles, 'roi': roi, 'ident': ident,
    'bands': bands_meta, 'per': per,
    'gel_png': 'data:image/png;base64,'+gel_b64,
}
with open('payload.json','w') as f:
    json.dump(payload, f, separators=(',',':'))
import os
print('payload bytes:', os.path.getsize('payload.json'))
print('gel png b64 KB:', round(len(gel_b64)/1024,1), 'crop', crop.shape, '-> disp', (target_w,target_h))
print('mm_per_px', round(mm_per_px,5), 'max mm', round((N-1)*mm_per_px,2))
print('lanes', lanes)
print('measured/ref lanes', sorted(int(k) for k in ident))
