function startSearch() {
    const headsup = document.getElementById("headsup");
    headsup.style.display = "none";

    const start = document.getElementById("start");
    start.style.display = "none";

    const inProgress = document.getElementById("inProgress");
    inProgress.style.display = "block";

    window.location.href = "/in_progress";
}
