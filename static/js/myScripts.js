function startSearch() {
    var headsup = document.getElementById("headsup");
    headsup.style.display = "none";

    var start = document.getElementById("start");
    start.style.display = "none";

    var inProgress = document.getElementById("inProgress");
    inProgress.style.display = "block";

    window.location.href = "/in_progress";
}
