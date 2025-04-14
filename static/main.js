
(function () {
  // Dark mode toggle
  const toggle = document.getElementById("darkModeToggle");
  if (toggle) {
    const darkPref = localStorage.getItem("dark-mode");
    if (darkPref === "true") {
      document.body.classList.add("dark-mode");
      toggle.checked = true;
    }

    toggle.addEventListener("change", () => {
      document.body.classList.toggle("dark-mode", toggle.checked);
      localStorage.setItem("dark-mode", toggle.checked ? "true" : "false");
    });
  }

  // Scroll to results if available
  window.addEventListener('DOMContentLoaded', function () {
    setTimeout(() => {
      const results = document.getElementById("results-section");
      if (results) {
        results.scrollIntoView({ behavior: "smooth" });
      }
    }, 300);

    // Show/hide back-to-top button
    const backToTop = document.getElementById("back-to-top");
    if (backToTop) {
      window.addEventListener('scroll', () => {
        backToTop.style.display = window.scrollY > 300 ? 'block' : 'none';
      });
    }
  });

  // Scroll to top
  window.scrollToTop = function () {
    window.scrollTo({ top: 0, behavior: 'smooth' });
  };

  // Toggle info box
  window.toggleInfo = function (id) {
    const el = document.getElementById(id);
    if (el) el.classList.toggle('info-text-show');
  };

  // Form validation
  const form = document.querySelector('form');
  if (form) {
    form.addEventListener("submit", function (event) {
      const seqInput = document.querySelector('textarea[name="sequence"]');
      const sequence = seqInput.value.trim().toUpperCase();
      const errorMessage = document.querySelector('#input-error');
      const sizeRangeInput = document.querySelector('input[name="product_size_range"]');
      let sizeRange = sizeRangeInput.value.trim() || '70-150';
      let minSize = 70;
      let rangeValid = true;

      if (errorMessage) errorMessage.remove();

      if (/^\d+-\d+$/.test(sizeRange)) {
        const [lower, upper] = sizeRange.split('-').map(Number);
        if (lower < 50 || upper <= lower) {
          rangeValid = false;
        } else {
          minSize = lower;
        }
      } else {
        rangeValid = false;
      }

      if (!rangeValid) {
        event.preventDefault();
        const msg = document.createElement("div");
        msg.id = "input-error";
        msg.className = "error-message";
        msg.textContent = "⚠ Please enter a valid product size range like '70-150' (min ≥ 50 and min < max).";
        sizeRangeInput.insertAdjacentElement("afterend", msg);
        return;
      }

      const isValid = /^[ACGT]+$/.test(sequence);
      if (!isValid || sequence.length < minSize) {
        event.preventDefault();
        const msg = document.createElement("div");
        msg.id = "input-error";
        msg.className = "error-message";
        msg.textContent = !isValid
          ? "⚠ Please enter a valid DNA sequence (A/T/C/G only)."
          : `⚠ Sequence must be at least ${minSize} bp to match the minimum product size.`;
        seqInput.insertAdjacentElement("afterend", msg);
      } else {
        seqInput.value = sequence; // Clean and uppercase
      }
    });
  }

  // Copy primer info
  window.copyPrimerBlock = function (button) {
    const container = button.closest(".primer-card");
    const content = container.querySelector(".primer-block-copy");
    if (content) {
      const text = content.textContent.trim();
      navigator.clipboard.writeText(text).then(() => {
        button.textContent = "✅";
        setTimeout(() => {
          button.textContent = "📋";
        }, 1500);
      }).catch(err => {
        console.error("Copy failed:", err);
        alert("Unable to copy. Please try manually.");
      });
    }
  };
})();
