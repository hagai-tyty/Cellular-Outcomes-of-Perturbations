from cellfate.common.progress import ProgressTracker


def test_mark_done_and_resume(tmp_path):
    p = tmp_path / "progress_tracker.json"
    t = ProgressTracker(p)
    assert not t.is_done("c1")
    t.mark_done("c1", 100)
    t.mark_done("c2", 50)
    assert t.is_done("c1") and t.n_done == 2 and t.n_samples == 150

    # a fresh tracker re-reads persisted state (simulates a Colab restart)
    t2 = ProgressTracker(p)
    assert t2.is_done("c1") and t2.is_done("c2")
    assert t2.pending(["c1", "c2", "c3"]) == ["c3"]


def test_failed_then_done_clears_failure(tmp_path):
    p = tmp_path / "progress_tracker.json"
    t = ProgressTracker(p)
    t.mark_failed("c1", "boom")
    assert "c1" in t.state.failed
    t.mark_done("c1", 10)
    assert "c1" not in t.state.failed and t.is_done("c1")
