/* Scroll effects for the landing page (overrides/home.html). Everything here is decoration: without JavaScript, or with
   "reduce motion" turned on, the page shows the same content as a static layout. */
(function () {
  var home = document.querySelector(".mf-home");
  if (!home) return;
  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) return;
  var root = document.documentElement;
  root.classList.add("mf-motion");

  var clamp = function (v) { return v < 0 ? 0 : v > 1 ? 1 : v; };
  var wide = window.matchMedia("(min-width: 60em)");

  /* Reading progress bar along the top of the window. */
  var bar = document.createElement("div");
  bar.className = "mf-progress";
  bar.setAttribute("aria-hidden", "true");
  document.body.appendChild(bar);

  /* Fade-up reveals, staggered within each group. */
  var groups = [
    [".mf-home h2"], [".mf-note"], [".mf-steps > li"], [".mf-check > li"], [".mf-files tr"],
    [".mf-pipe > li"], [".mf-caps > div"], [".mf-lipids > figure"], [".mf-home details.mf-more"]
  ];
  var revealed = [];
  groups.forEach(function (g) {
    home.querySelectorAll(g[0]).forEach(function (el, i) {
      el.classList.add("mf-reveal");
      el.style.setProperty("--mf-delay", Math.min(i, 5) * 70 + "ms");
      revealed.push(el);
    });
  });
  if ("IntersectionObserver" in window) {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) { e.target.classList.add("is-in"); io.unobserve(e.target); }
      });
    }, { rootMargin: "0px 0px -8% 0px", threshold: 0.08 });
    revealed.forEach(function (el) { io.observe(el); });
  } else {
    revealed.forEach(function (el) { el.classList.add("is-in"); });
  }

  /* Hero: the membrane grows toward the reader as the hero scrolls away. */
  var hero = home.querySelector(".mf-hero");
  var heroMedia = home.querySelector(".mf-hero-media");

  /* "From structure to simulation": on wide screens the four stages become one pinned panel. Scrolling steps through
     them, and each picture grows a little while it is on screen. On narrow screens the grid stays and each picture
     zooms in as it enters. */
  var stages = home.querySelector("[data-mf-story]");
  var figs = stages ? Array.prototype.slice.call(stages.querySelectorAll("figure")) : [];
  var pin = null, dots = [];
  function buildStory() {
    if (!stages || pin) return;
    pin = document.createElement("div");
    pin.className = "mf-story-pin";
    figs.forEach(function (f) { pin.appendChild(f); });
    var nav = document.createElement("ol");
    nav.className = "mf-story-dots";
    nav.setAttribute("aria-hidden", "true");
    figs.forEach(function (f, i) {
      var li = document.createElement("li");
      var b = f.querySelector("figcaption b");
      li.textContent = b ? b.textContent.replace(/^\d+\.\s*/, "").replace(/\.$/, "") : "Step " + (i + 1);
      nav.appendChild(li);
      dots.push(li);
    });
    pin.appendChild(nav);
    stages.appendChild(pin);
    stages.classList.add("mf-story");
    stages.style.setProperty("--mf-n", figs.length);
  }
  function dropStory() {
    if (!pin) return;
    figs.forEach(function (f) { stages.insertBefore(f, pin); f.classList.remove("is-active"); f.style.removeProperty("--mf-z"); });
    pin.remove(); pin = null; dots = [];
    stages.classList.remove("mf-story");
  }
  function layout() { if (wide.matches) buildStory(); else dropStory(); update(); }
  if (wide.addEventListener) wide.addEventListener("change", layout); else wide.addListener(layout);

  /* MD movie: frames are fetched when the section comes near, then scrolling picks the frame and zooms in. */
  var movie = home.querySelector("[data-mf-movie]");
  var movieImg = movie ? movie.querySelector(".mf-movie-frame img") : null;
  var movieFrames = [], movieShown = -1;
  if (movie && movieImg) {
    var dir = movie.getAttribute("data-frame-dir"), count = +movie.getAttribute("data-frame-count") || 0;
    var load = function () {
      movieFrames = []; movieShown = -1;
      for (var i = 0; i < count; i++) {
        var im = new Image();
        im.src = dir + "frame_" + (i < 10 ? "0" : "") + i + ".jpg";
        movieFrames.push(im);
      }
    };
    if ("IntersectionObserver" in window) {
      var mo = new IntersectionObserver(function (es) {
        if (es[0].isIntersecting) { load(); update(); mo.disconnect(); }
      }, { rootMargin: "150% 0px" });
      mo.observe(movie);
    } else { load(); }
  }

  var ticking = false;
  function update() {
    ticking = false;
    var vh = window.innerHeight;
    var max = root.scrollHeight - vh;
    bar.style.transform = "scaleX(" + (max > 0 ? clamp(window.scrollY / max) : 0) + ")";

    if (hero && heroMedia) {
      var r = hero.getBoundingClientRect();
      var p = clamp(-r.top / Math.max(r.height, 1));
      heroMedia.style.setProperty("--mf-hp", p.toFixed(4));
    }

    if (pin) {
      var s = stages.getBoundingClientRect();
      var run = s.height - pin.offsetHeight;
      var prog = clamp(-s.top / Math.max(run, 1));
      var n = figs.length;
      var idx = Math.min(n - 1, Math.floor(prog * n));
      figs.forEach(function (f, i) {
        var local = clamp(prog * n - i);
        f.classList.toggle("is-active", i === idx);
        f.style.setProperty("--mf-z", (0.9 + 0.14 * local).toFixed(4));
      });
      dots.forEach(function (d, i) { d.classList.toggle("is-active", i === idx); d.classList.toggle("is-done", i < idx); });
    } else if (figs.length) {
      figs.forEach(function (f) {
        var b = f.getBoundingClientRect();
        var t = clamp((vh - b.top) / (vh * 0.7));
        f.style.setProperty("--mf-z", (0.86 + 0.14 * t).toFixed(4));
      });
    }

    if (movie && movieImg) {
      var m = movie.getBoundingClientRect();
      var mp = clamp(-m.top / Math.max(m.height - vh, 1));
      movie.style.setProperty("--mf-mp", mp.toFixed(4));
      movie.style.setProperty("--mf-mz", (0.82 + 0.38 * mp).toFixed(4));
      if (movieFrames.length) {
        var k = Math.min(movieFrames.length - 1, Math.round(mp * (movieFrames.length - 1)));
        if (k !== movieShown && movieFrames[k].complete && movieFrames[k].naturalWidth) { movieImg.src = movieFrames[k].src; movieShown = k; }
      }
    }
  }
  function onScroll() { if (!ticking) { ticking = true; window.requestAnimationFrame(update); } }
  window.addEventListener("scroll", onScroll, { passive: true });
  window.addEventListener("resize", onScroll);
  layout();
})();
