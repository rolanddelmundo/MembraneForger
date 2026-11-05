/* Landing page behaviour (overrides/home.html): full-screen view of any render, the turning movie, a reading progress
   bar and fade-up reveals. Without JavaScript the page shows the same pictures as a static layout; with "reduce motion"
   on, only the full-screen view is kept. */
(function () {
  var home = document.querySelector(".mf-home");
  if (!home) return;
  var root = document.documentElement;
  var clamp = function (v) { return v < 0 ? 0 : v > 1 ? 1 : v; };

  /* A picture that does not exist yet (a lipid render not committed) hides its figure instead of showing a broken image. */
  home.querySelectorAll(".mf-lipids img").forEach(function (im) {
    var hide = function () { var f = im.closest("figure"); if (f) f.style.display = "none"; };
    if (im.complete && im.naturalWidth === 0 && im.getAttribute("src")) hide();
    im.addEventListener("error", hide);
  });

  /* Full-screen view: select any render to open its largest copy. */
  var box = home.querySelector(".mf-lightbox");
  if (box && box.showModal) {
    var boxImg = box.querySelector("img"), boxCap = box.querySelector(".mf-lightbox-cap");
    home.querySelectorAll(".mf-zoom img").forEach(function (im) {
      im.setAttribute("tabindex", "0");
      im.setAttribute("role", "button");
      var open = function () {
        var fig = im.closest("figure");
        var cap = fig ? fig.querySelector("figcaption") : null;
        boxImg.src = im.currentSrc || im.src;
        var full = im.getAttribute("data-full");
        if (full) { var pre = new Image(); pre.onload = function () { boxImg.src = full; }; pre.src = full; }
        boxImg.alt = im.alt;
        boxCap.textContent = cap ? cap.textContent.trim() : im.alt;
        box.showModal();
      };
      im.addEventListener("click", open);
      im.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); open(); } });
    });
    box.addEventListener("click", function () { box.close(); });
  }

  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) return;
  root.classList.add("mf-motion");

  var bar = document.createElement("div");
  bar.className = "mf-progress";
  bar.setAttribute("aria-hidden", "true");
  document.body.appendChild(bar);

  var revealed = [];
  [".mf-home h2", ".mf-note", ".mf-strip > figure", ".mf-pipe > li", ".mf-cards > div", ".mf-lipids > figure",
   ".mf-stats > div", ".mf-reqs > a", ".mf-next > a", ".mf-check > li"].forEach(function (sel) {
    home.querySelectorAll(sel).forEach(function (el, i) {
      el.classList.add("mf-reveal");
      el.style.setProperty("--mf-delay", Math.min(i, 6) * 60 + "ms");
      revealed.push(el);
    });
  });
  if ("IntersectionObserver" in window) {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { if (e.isIntersecting) { e.target.classList.add("is-in"); io.unobserve(e.target); } });
    }, { rootMargin: "0px 0px -6% 0px", threshold: 0.05 });
    revealed.forEach(function (el) { io.observe(el); });
  } else {
    revealed.forEach(function (el) { el.classList.add("is-in"); });
  }

  /* Movie: scrolling picks the frame. The frames are fetched when the section comes near, from the first directory in
     data-frame-dirs that has them (the 1920 px copies first). */
  var movie = home.querySelector("[data-mf-movie]");
  var movieImg = movie ? movie.querySelector(".mf-movie-frame img") : null;
  var frames = [], shown = -1;
  if (movie && movieImg) {
    var dirs = (movie.getAttribute("data-frame-dirs") || "").split(/\s+/).filter(Boolean);
    var count = +movie.getAttribute("data-frame-count") || 0;
    var name = function (i) { return "frame_" + (i < 10 ? "0" : "") + i + ".jpg"; };
    var tryDir = function (k) {
      if (k >= dirs.length) return;
      var probe = new Image();
      probe.onload = function () {
        movieImg.src = probe.src;
        frames = []; shown = 0;
        for (var i = 0; i < count; i++) { var im = new Image(); im.src = dirs[k] + name(i); frames.push(im); }
        update();
      };
      probe.onerror = function () { tryDir(k + 1); };
      probe.src = dirs[k] + name(0);
    };
    if ("IntersectionObserver" in window) {
      var mo = new IntersectionObserver(function (es) { if (es[0].isIntersecting) { tryDir(0); mo.disconnect(); } }, { rootMargin: "150% 0px" });
      mo.observe(movie);
    } else { tryDir(0); }
  }

  var ticking = false;
  function update() {
    ticking = false;
    var vh = window.innerHeight, max = root.scrollHeight - vh;
    bar.style.transform = "scaleX(" + (max > 0 ? clamp(window.scrollY / max) : 0) + ")";
    if (movie && movieImg) {
      var m = movie.getBoundingClientRect();
      var mp = clamp(-m.top / Math.max(m.height - vh, 1));
      movie.style.setProperty("--mf-mp", mp.toFixed(4));
      if (frames.length) {
        var k = Math.min(frames.length - 1, Math.round(mp * (frames.length - 1)));
        if (k !== shown && frames[k].complete && frames[k].naturalWidth) { movieImg.src = frames[k].src; shown = k; }
      }
    }
  }
  function onScroll() { if (!ticking) { ticking = true; window.requestAnimationFrame(update); } }
  window.addEventListener("scroll", onScroll, { passive: true });
  window.addEventListener("resize", onScroll);
  update();
})();
