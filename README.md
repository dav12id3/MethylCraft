# MethylCraft

##  Demo Video
[Watch the demo on YouTube](https://youtu.be/jaD5FA22yLk)

##  Project Description

**MethylCraft** is a responsive, interactive web tool for designing methylation PCR primers and methylation-specific probes for bisulfite-converted DNA sequences. It is desgined for researchers and students investigating DNA methylation.

Built on top of the `primer3-py` library (a Python API for Primer3), the app allows users to:

- Automatically perform **in silico bisulfite conversion**
- Customize primer/probe design: Tm, GC content, product size, salt/dNTP concentrations
- Identify both methylation-sensitive (M-Probe) and methylation-insensitive (U-Probe) probes
- Quickly copy sequences and view detailed thermodynamic metrics

## Detailed desrciption for project structure and individual files

Below is a description of the key files and folders that make up the MethylCraft web application:

```bash
MethylCraft/
├── app.py
├── helper.py
├── templates/
│   └── index.html
├── static/
│   ├── style.css
│   └── main.js
├── README.md
├── requirements.txt
```

---

### 🔧 Backend

| File         | Description |
|--------------|-------------|
| **`app.py`** | The main Flask application. Handles routing (`GET`) when users visit websiite and (`POST`) when users input sequence and searching criteria. It also perform input validation, in silico bisulfite conversion, primer/probe generation via `primer3-py`, and rendering of results using `index.html`. |
| **`helper.py`** | Contains core helper functions used by `app.py`, including: <br>– `bisulfite_convert()` for DNA transformation <br>– `run_primer3()` for calling Primer3 using user inputs (if provided) <br>– `design_probe_from_fixed_primers()` for U-Probe finding <br>– Other functions such as CpG counting, degenerate primer creation, probe filtering, and highlighting logic are also defined and called in the core hekper functions. |

---

### 🎨 Frontend

| File          | Description |
|---------------|-------------|
| **`templates/index.html`** | The main HTML page rendered by Flask. Includes form inputs, results display (primer/probe cards), interactive tooltips, and a dynamic UI. Uses Jinja2 templating for inserting results. |
| **`static/style.css`**     | Custom CSS for layout, dark mode, responsive design, tooltips, and primer/probe card formatting. Includes sticky footer, form alignment, and themed UI styling. |
| **`static/main.js`**       | JavaScript for client-side features: <br>– Dark mode toggle <br>– Scroll-to-results and scroll-to-top <br>– Form input validation <br>– Copy-to-clipboard <br>– Tooltip toggle behavior. |

---

### 📝 Other Files

| File             | Description |
|------------------|-------------|
| **`requirements.txt`** | Lists Python dependencies required by Render or any Python environment (e.g., Flask, primer3-py, gunicorn if deployed). |