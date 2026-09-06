set -euo pipefail
python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
old='''            items.append(
                f'{row["donor_id"]} (err={row["error"]:+.4f}, '
                f'braak={int(row["braak_numeric"]) if pd.notna(row["braak_numeric"]) else "NA"}, '
                f'cerad={int(row["cerad_ordinal"]) if pd.notna(row["cerad_ordinal"]) else "NA"}, '
                f'{row["cognitive_status"]})'
            )'''
new='''            items.append(
                f'{row["donor_id"]} (err={row["error"]:+.4f}, '
                f'braak={int(row["braak_numeric"]) if ("braak_numeric" in row.index and pd.notna(row["braak_numeric"])) else "NA"}, '
                f'cerad={int(row["cerad_ordinal"]) if ("cerad_ordinal" in row.index and pd.notna(row["cerad_ordinal"])) else "NA"}, '
                f'sex={row["sex"] if "sex" in row.index else "NA"})'
            )'''
text=text.replace(old,new)
path.write_text(text)
print('patched donor_list')
PY
python /scratch/result.py
echo
echo '--- results.json ---'
sed -n '1,220p' /scratch/results.json
echo
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
