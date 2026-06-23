from spectral_ibr.infrastructure.andes.state_partitioner import StatePartitioner


def test_state_partitioner_separates_control_markers() -> None:
    partition = StatePartitioner().partition(
        ("Bus34.delta", "PLL1.theta", "REGF1.xp", "Line12.current")
    )

    assert partition.electrical == ("Bus34.delta", "Line12.current")
    assert partition.control == ("PLL1.theta", "REGF1.xp")
