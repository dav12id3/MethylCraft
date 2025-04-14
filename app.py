from flask import Flask, render_template, request, flash, redirect, url_for
from helper import bisulfite_convert, design_primers # type: ignore

app = Flask(__name__)
app.secret_key = 'my_super_secret_key_123'

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':

        sequence = request.form.get("sequence", "").strip().upper()

        # Product size range
        size_range = request.form.get('product_size_range') or '70-150'
        # Check product size input format
        try:
            min_size, max_size = map(int, size_range.split('-'))
            product_range = [[min_size, max_size]]
        except:
            error_message = "❌ Invalid product size range format. Please use something like '70-150'."
            return render_template('index.html', input=sequence, error=error_message)

        if len(sequence) < min_size:
            flash(f"Sequence is too short (minimum required: {min_size} bp).")
            return redirect(url_for('index'))


        # Optional parameters from user input
        primer_min_tm = request.form.get('primer_min_tm') or 52.0
        primer_opt_tm = request.form.get('primer_opt_tm') or 55.0
        primer_max_tm = request.form.get('primer_max_tm') or 58.0

        probe_min_tm = request.form.get('probe_min_tm') or 57.0
        probe_opt_tm = request.form.get('probe_opt_tm') or 60.0
        probe_max_tm = request.form.get('probe_max_tm') or 63.0

        salt_mono = request.form.get('salt_mono') or 50.0
        salt_div = request.form.get('salt_div') or 3.0
        dntp_conc = request.form.get('dntp_conc') or 0.8

        # Perform in silico bisulfite conversion
        methylated_seq, unmethylated_seq = bisulfite_convert(sequence) 

        # Package parameters into a config dictionary
        primer_config = {
            'primer_min_tm': float(primer_min_tm),
            'primer_opt_tm': float(primer_opt_tm),
            'primer_max_tm': float(primer_max_tm),
            'probe_min_tm': float(probe_min_tm),
            'probe_opt_tm': float(probe_opt_tm),
            'probe_max_tm': float(probe_max_tm),
            'salt_mono': float(salt_mono),
            'salt_div': float(salt_div),
            'dntp_conc': float(dntp_conc),
            'product_size_range': product_range
        }

        methylated_primers = design_primers(
            methylated_seq,
            unmethylated_seq=unmethylated_seq,
            original_seq=sequence,
            config=primer_config
        )

        return render_template('index.html',
                               result={'methylated_primers': methylated_primers},
                               input=sequence,
                               submitted=True)

    return render_template('index.html', result=None, submitted=False)


if __name__ == '__main__':
    app.run(debug=True)
